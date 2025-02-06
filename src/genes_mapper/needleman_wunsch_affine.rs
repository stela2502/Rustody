const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -2;
const GAP_OPEN_PENALTY: i32 = -5; // Penalty for opening a gap
const GAP_EXTENSION_PENALTY: i32 = -4; // Penalty for extending a gap

use crate::genes_mapper::cigar::{Cigar, CigarEnum, CigarTuple};
use crate::traits::BinaryMatcher;

use std::fs::File;
use std::io::{self,Write};


pub struct NeedlemanWunschAffine{
	dp: Vec<Vec<i32>>,
	n: usize,
	m: usize,
	pub cigar: Cigar,
	//circles:usize,
	debug:bool,
}



impl <'a> NeedlemanWunschAffine {
	pub fn new( ) ->Self {
		// Initialize the DP matrix

		let me = Self {
	    	dp : Vec::<Vec::<i32>>::new(),
	    	n: 0,
	    	m: 0,
	    	cigar: Cigar::default(),
	    	//circles:0,
	    	debug: false,
	    };
	    me
	}

	pub fn cigar(&self) -> Cigar{
		self.cigar.clone()
	}

	pub fn debug( &self ) -> bool{
		self.debug
	}
	pub fn set_debug(&mut self, state:bool){
		self.debug = state
	}

	/// initialize does exactly that - initializes the internal storage - if necessary.
	pub fn initialize( &mut self, rows:usize, cols:usize ){

		let size = rows.max(cols);
		self.cigar.clear();

		if self.dp.len() < size +2 || self.dp[0].len() < size +2 {
			self.dp = vec![vec![0; size + 2]; size + 2];
			for i in 1..=size {
		        self.dp[i][0] = self.dp[i - 1][0] + GAP_OPEN_PENALTY + (i - 1) as i32 * GAP_EXTENSION_PENALTY;
		        self.dp[0][i] = self.dp[0][i - 1] + GAP_OPEN_PENALTY + (i - 1) as i32 * GAP_EXTENSION_PENALTY;
		    }
		}

	}

	pub fn to_string<T>( &mut self, read: &T, database: &T, _humming_cut: f32 ) -> String 
	where
    T: BinaryMatcher + std::fmt::Display{
    	format!("alignement:\n{}",self.cigar.as_alignement(read, database ))
	}

	pub fn int_state_to_string<T>( &mut self, read: &T, database: &T, cigar_vec: &[CigarTuple] ) -> String 
	where
    T: BinaryMatcher + std::fmt::Display{
    	self.cigar.reset_fom_path( cigar_vec );

    	format!("alignement: {}", self.cigar.as_alignement( read, database))
	}

	/// calculates a needleman_wunsch_affine matrix and returns the final value scaled to the shorter sequence length.
	pub fn needleman_wunsch_affine<T>(&mut self, read: &T, database: &T, humming_cut: f32 ) -> f32 
	where
    T: BinaryMatcher {

    	if read.as_dna_string() == database.as_dna_string() {
    		// a 100% match?!
    		self.cigar.restart_from_cigar( &format!("{}M", database.len()) );
    		return 0.0
    	}

	    self.n = read.len();
	    self.m = database.len();
	    let n = self.n;
	    let m = self.m;
	    if n.min(m) < 15  && read.tri_nuc_abs_diff(database) > humming_cut{
            return 100.0
        }

        self.initialize( n, m );

		// Fill in the DP matrix
		for i in 1..=n {
		    for j in 1..=m {
		        if let (Some(nuc1), Some(nuc2)) = (read.get_nucleotide_2bit(i - 1), database.get_nucleotide_2bit(j - 1)) {
		            // Calculate match/mismatch score
		            let match_mismatch_score = if nuc1 == nuc2 { MATCH_SCORE } else { MISMATCH_SCORE };

		            let gap_open_penalty = if match_mismatch_score == MISMATCH_SCORE && read.is_same_streak(i - 1) {
		            	MATCH_SCORE
		            }else {
						GAP_OPEN_PENALTY
		            };
		            //let gap_open_penalty = GAP_OPEN_PENALTY;

		            // Compute the maximum score directly
		            self.dp[i][j] = [
		            	// Score for diagonal (match/mismatch)
		                self.dp[i - 1][j - 1] + match_mismatch_score,
		                // Score for gap in read (insertion in read)
		                self.dp[i - 1][j] + gap_open_penalty,
		                // Score for gap in database (deletion in read)
		                self.dp[i][j - 1] + gap_open_penalty, 
		                // Score for gap extension in read (insertion in read)
		                (1..=i - 1).map(|k| self.dp[i - k][j] + gap_open_penalty + k as i32 * GAP_EXTENSION_PENALTY).max().unwrap_or(std::i32::MIN),
		                // Score for gap extension in database (deletion in read)
		                (1..=j - 1).map(|k| self.dp[i][j - k] + gap_open_penalty + k as i32 * GAP_EXTENSION_PENALTY).max().unwrap_or(std::i32::MIN),
		            ].iter().copied().max().unwrap();
		        } else {
		            // Panic if unable to retrieve nucleotide sequences
		            panic!("Could not get sequence for nuc1 or nuc2 needleman_wunsch_affine")
		        }
		    }
		}
	    let size = n.min(m) as f32;
	    self.to_cigar( read, database, humming_cut );
	    // this call might have swapped but the self.n and self.m have been updated!
	    (size - self.dp[self.n][self.m] as f32) / size

	}

	/// This will give two vectors - the aligned sequence 1 and sequence 2
	pub fn needleman_wunsch_affine_backtrack<T>(&mut self, read: &T, database: &T, cig_vec: &[CigarEnum] ) -> ( String, String ) 
	where
    T: BinaryMatcher {
		//let (mut i, mut j) = (read.len(), database.len());


		let (mut i, mut j);
	    let mut align_read = Vec::<u8>::with_capacity( cig_vec.len() );
	    let mut align_database = Vec::<u8>::with_capacity( cig_vec.len() );
	    i = 0;
	    j = 0;

	    for val in cig_vec {
	    	match val {
	    		CigarEnum::Insertion => {
	    			if let Some(nuc1) = read.get_nucleotide_2bit(i){
	    				align_read.push( Self::to_utf8_lc(nuc1) );
	            		align_database.push( b'-' );
	            		i += 1;
	    			}else {
	    				panic!("Could not get sequence for nuc1 #{i}")
	    			}
	    		},
	    		CigarEnum::Deletion => {
	    			if let Some(nuc2) = database.get_nucleotide_2bit(j){
	    				align_read.push( b'-' ); 
						align_database.push( Self::to_utf8_lc(nuc2) );
						j +=1;
	    			}else {
	    				eprintln!("Could not get sequence for nuc2 #{j} Something is seriousely wrong here!!!");
	    				break;
	    			}
	            	
	    		},
	    		CigarEnum::Mismatch => {
	    			if let Some(nuc1) =  read.get_nucleotide_2bit(i){
	    				align_read.push( Self::to_utf8_lc( nuc1 ) );
	    				i += 1;
	    			}else {
	    				align_read.push( b'-' ); 
	    			}

	    			if let Some(nuc2) =  database.get_nucleotide_2bit(j){
	    				align_database.push( Self::to_utf8_lc( nuc2 ) );
	    				j += 1;
	    			}else {
	    				align_database.push( b'-' ); 
	    			}
	    		},
	    		_ => {

	    			if let Some(nuc1) =  read.get_nucleotide_2bit(i){
	    				align_read.push( Self::to_utf8( nuc1 ) );
	    				i += 1;
	    			}else {
	    				align_read.push( b'-' ); 
	    			}
	    			if let Some(nuc2) =  database.get_nucleotide_2bit(j){
	    				align_database.push( Self::to_utf8( nuc2 ) );
	    				j += 1;
	    			}else {
	    				align_database.push( b'-' ); 
	    			}
	    		},
	    	}

	    }

	    ( String::from_utf8_lossy( &align_read).to_string(), String::from_utf8_lossy( &align_database).to_string() )
	}

	fn to_utf8( twobit:u8 ) -> u8 {
		match twobit{
			0 => b'A',
			1 => b'C',
			2 => b'G',
			3 => b'T',
			_ => b'N',
		}
	}

	fn to_utf8_lc( twobit:u8 ) -> u8 {
		match twobit{
			0 => b'a',
			1 => b'c',
			2 => b'g',
			3 => b't',
			_ => b'n',
		}
	}

	/*pub fn to_cigar_vec_inverted<T>( &mut self, read: &T, database: &T, humming_cut: f32 )-> Vec<CigarEnum> 
	where
    T: BinaryMatcher {

    	self.circles +=1;
    	if  self.circles < 2 {
    		let _ =self.needleman_wunsch_affine( database, read, humming_cut );
    	}
	    let results = self.cigar_vec();
	    self.circles = 0;

	    results.iter()
    		.map(|cigar_enum| if *cigar_enum == CigarEnum::Deletion { CigarEnum::Insertion } else if *cigar_enum == CigarEnum::Insertion { CigarEnum::Deletion } else { *cigar_enum })
    		.collect()
	}*/

	pub fn cigar_vec( &self ) -> Vec<CigarEnum>{
		self.cigar.to_vec()
	}


	/// This will create the Cigar states vector from the internal matrix.
	pub fn to_cigar<T>(&mut self, read: &T, database: &T, _humming_cut: f32 ) 
	where
    T: BinaryMatcher {

    	//println!("needleman_wunsch_affine::to_cigar_vec is called.");

    	//self.cigar.calculate_cigar( &self.dp, read.get_nucleotide_2bit(read.len()-1) == database.get_nucleotide_2bit(database.len() -1) );
		
		let (mut i, mut j) = (read.len(), database.len());

		if self.cigar.len() == i.max(j) {
			// this seams to have been run already?!
			return;
		}
			
		let mut path = Vec::<CigarTuple>::with_capacity(i.max(j));	    
		//let mut rev_id = i.max(j).saturating_sub(1);

		//if read.get_nucleotide_2bit(0) == database.get_nucleotide_2bit(0)

	    while i > 0 || j > 0 {

	    	//if let (Some(nuc1), Some(nuc2)) = (read.get_nucleotide_2bit(i.saturating_sub(1)), database.get_nucleotide_2bit(j.saturating_sub(1)) ){
    		let this_value = self.dp[i][j];
    		#[cfg( debug_assertions )]
    		if self.debug {
    			println!("I am testing this part of the dp matrix with the lower right corner i:{i};j:{j}:\n{}\t{}\n{}\t{}",
	    			self.dp[i.saturating_sub(1)][j.saturating_sub(1)], self.dp[i][j.saturating_sub(1)],
	    			self.dp[i.saturating_sub(1)][j],self.dp[i][j]
	    			);
    		}
	    	if let Some((max_index, max_value)) = vec![
	    	self.dp[i.saturating_sub(1)][j.saturating_sub(1)], 
	    	self.dp[i.saturating_sub(1)][j],
	    	self.dp[i][j.saturating_sub(1)]
	    	].iter().enumerate().max_by_key(|(_, &val)| val) {
	    		let option = match max_index{
	    		//cigar[rev_id] = match max_index{
	    			0 => {
	    				// a match or a mismatch self.dp[i.saturating_sub(1)][j.saturating_sub(1)]
	    				i = i.saturating_sub(1);
	    				j =j.saturating_sub(1);
	    				if  *max_value == this_value - MATCH_SCORE {
	    					CigarEnum::Match
	    				}else {
	    					CigarEnum::Mismatch
	    				}
	    			},
	    			1=>{
	    				// an insertion is most likely
	    				i = i.saturating_sub( 1 );
	    				CigarEnum::Insertion
	    			},
	    			2=> {
	    				// a deletion is most likely
	    				j = j.saturating_sub(1);
	    				CigarEnum::Deletion
	    			},
	    			_ => unreachable!()
	    		};
	    		if let Some(tuple) = path.last_mut() {
				    if tuple.is_a(&option) {
				        tuple.vec_len += 1;
				    } else {
				    	//println!("After {} we start with an option {}", tuple, option);
				        path.push(CigarTuple::from_scratch(option, 1));
				    }
				} else {
				    path.push(CigarTuple::from_scratch(option, 1));
				}
	    		#[cfg( debug_assertions )]
	    		if self.debug{
	    			println!("inserted one additional to the last elemet {:?}",  path.last() );
	    		}
	    		/*println!("Here ({i};{j} I had a max index of {max_index} and a max value of {max_value} and decided on a {} with read:{:?} and query:{:?}",
	    			cigar[rev_id], read.get_nucleotide_2bit(i), database.get_nucleotide_2bit(j) );
	    		println!("This is the the matrix that lead to that:\n\n{}\t{}\n{}\t*{}*",
	    			self.dp[i.saturating_sub(1)][j.saturating_sub(1)], self.dp[i][j.saturating_sub(1)],
	    			self.dp[i.saturating_sub(1)][j],self.dp[i][j] );
	    		*/
	    		
	    	}else {
	    		panic!("I can not decode the df matrix to a Cigar state at :{i}; j{j}:\n{}\t{}\n{}\t*{}*",
	    			self.dp[i.saturating_sub(1)][j.saturating_sub(1)], self.dp[i][j.saturating_sub(1)],
	    			self.dp[i.saturating_sub(1)][j],self.dp[i][j]
	    			);
	    	}
		    /*if rev_id > 0 {
		    	rev_id -= 1;
		    }else if rev_id == 0{
		    	if i > 0 {
			    	let option = CigarEnum::Insertion;
			    	if let Some(tuple) = path.last_mut() {
					    if tuple.is_a(&option) {
					        tuple.vec_len += i;
					    } else {
					        path.push(CigarTuple::from_scratch(option, 1));
					    }
					} else {
					    unreachable!();
					}
				}if j > 0 {
			    	//println!("End not reached by main - adding DELETION");
			    	let option = CigarEnum::Deletion;
			    	if let Some(tuple) = path.last_mut() {
					    if tuple.is_a(&option) {
					        tuple.vec_len += j;
					    } else {
					        path.push(CigarTuple::from_scratch(option, 1));
					    }
					} else {
					    unreachable!();
					}
			    }
			    break;
			} */
	    }

	    //if rev_id != 0 {
	    //	panic!("We have not filled in all the values here!?!?");
	    //}
	    let rev: Vec<_> = path.into_iter().rev().collect();
		self.cigar.reset_fom_path ( &rev );
		
		#[cfg(all(debug_assertions, feature = "mapping_debug"))]
		{
			println!("        nascent cigar string:\n{}", self.cigar.as_alignement(read, database ));
	    }
		self.cigar.fix_di_problems (0, read, database );

		self.cigar.finalize();
		self.cigar.check_alignment( read, database);
		( self.cigar.dropped_start, self.cigar.dropped_end ) = read.get_dropped_values();
		
		
		#[cfg(all(debug_assertions, feature = "mapping_debug"))]
		{	
			let mut cig = self.cigar.clone();
			println!("del/ins remapped cigar string:\n{}\nstate changes: {}; mapping quality {}\nSam Cigar: {}", 
				self.cigar.as_alignement(read, database ), 
				self.cigar.state_changes(),
				self.cigar.mapping_quality(),
				cig.to_sam_string().0,
				);
	    }

	}

	/// a rather useless function converting the Cigar vector to a String - No Cigar string - but the exptended version.
	pub fn cigar_to_string(cigar:&[CigarEnum]) -> String{
		let mut ret = "".to_string();
		for i in cigar{
			ret += &format!("{}",i)
		}
		ret
	}



	// Function to export DP matrix to a file
	pub fn export_dp_matrix(&self, file_path: &str) -> io::Result<()> {
	    let mut file = match File::create(file_path) {
	    	Ok(file) => file,
	    	Err(err) => {
	    		eprintln!( "Trying to write th dp matrix of a needleman wunsch affine mapping I hit this file system error:\n{err:?}");
	    		// this is not too crucial here and it would be a shame to kill the whole process over this!
	    		return Ok(()) 
	    	},
	    };

	    for row in &self.dp {
	        let row_str: Vec<String> = row.iter().map(|&val| val.to_string()).collect();
	        let row_line = row_str.join("\t");
	        writeln!(file, "{}", row_line)?;
	    }
	    println!("I have exported the needleman wunsch affine primary matrix to {file_path}");
	    Ok(())
	}


}


