//gene_mapper::cigar.rs

use crate::traits::Cell;
use crate::traits::BinaryMatcher;
//use crate::genes_mapper::gene_data::GeneData;
use crate::genes_mapper::cigar::{ CigarTuple, CigarEnum };


//use crate::genes_mapper::gene_data::GeneData;
use core::cmp::Ordering;


use regex::Regex;

use core::cmp::max;
use core::fmt;

//use std::fs::File;
//use std::io::Write;



/// This enum denotes the way a Cigar's end areas have been modified:
/// NA - no modification; Start - a longer caotic area at the start has been fixed -> end match;
/// End - a longer fix at the end of the Cigar has been fixed - a start match;
/// Both - and internal match?! Let's see if that even happens.
/// StartInsert - this cigar indicated that teh match was longer on the starting end - Likely the db too short?
///               Anyhow - this need to be fixed in the analysis scripts and therefore a CigarEndFix is needed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarEndFix{
	Na,
	Start,
	End,
	Both,
	StartInsert,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cigar{
	// the cigar string
	pub cigar:String,
	pub fixed:Option<CigarEndFix>,
	debug:bool,
	/// contains whether <CigarEnum>.to_id() is part of this cigar
	contains: Vec<bool>,
	/// state_changes is a measure of quality for a Cigar - the more changes the less likely that this is real.
	state_changes: usize,
	/// keep track of how many bp have been sliced from this entry's start
	pub dropped_start:usize,
	/// keep track of how many bp have been sliced from this entry's end
	pub dropped_end:usize,
	/// should be another test - max match > (self.len() as f16 * 0.7) as usize
	pub max_match: usize,

}

// Implementing Display trait for SecondSeq
impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Cigar {} - fixed {:?}; state changed {}; contains {:?} slice info {:?}", self.cigar, self.fixed, self.state_changes, self.contains ,( self.dropped_start, self.dropped_end) )
    }
}

impl Default for Cigar {
    fn default() -> Self {
        Cigar {
            cigar: "".to_string(),
            fixed: None,
            debug: false,
            contains: vec![false;7],
            state_changes: 1000,
            dropped_start: 0,
            dropped_end:0,
            max_match: 0,     
        }
    }
}

impl PartialOrd for Cigar {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cigar {
    fn cmp(&self, other: &Self) -> Ordering {
        // Define your comparison logic here. For example:
        // 1. Compare by `state_changes` (ascending)
        // 2. Then by length of `cigar` string (descending)
        // 3. Finally by `cigar` string lexicographically (ascending)

        self.mapping_quality().cmp(&other.mapping_quality()) // more is better
        .then_with( || self.len().cmp(&other.len())) // more is better
        .then_with( || other.state_changes.cmp(&self.state_changes)) // less is better
    }
}

impl Cigar{

	pub fn new( cigar: &str ) -> Self{

		let mut ret = Self::default();
		ret.cigar = cigar.to_string();
		let vec = ret.to_cigar_tupel_vec(false);

		ret.reset_fom_path( &vec );
		ret
	}

		/// converts a CigarEnum vector into a Cigar string and stores that internally.
	/// This function also updated the contains vector.
	pub fn reset_fom_path(&mut self, path: &[CigarTuple] ){

		self.cigar.clear();

		// boost and compact again:
		//let corrected = Self::finalize( path );
		self.state_changes = path.len();

	    for tuple in path {
	        self.cigar.push_str(&tuple.to_string());
	        self.contains[tuple.option.to_id()] = true;
	        if tuple.len() > self.max_match{
				self.max_match = tuple.len();
			}
	    }
	}

	pub fn is_decent ( &self ) -> bool {
		self.mapping_quality() > 30
		&& self.state_changes < 5
		|| self.max_match > (self.len() as f32 * 0.7) as usize
		
	}

	pub fn better_as( &self, other: &Self ) -> bool{
		self > other
	}

	pub fn set_debug(&mut self, state:bool) {
		self.debug = state;
	}

	pub fn state_changes(&self) -> usize{
		self.state_changes
	}
	
	pub fn clear(&mut self) {
		self.cigar ="".to_string();
		self.fixed = None;
		self.contains.fill(false);
		self.state_changes = 0;
		//self.dropped_start = 0;
		//self.dropped_end = 0;
	}

	fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

    pub fn to_string(&self) -> String {
    	self.cigar.to_string()
    }

    fn get_nucleotide_2bit( &self, bit: Option<u8> ) -> String {
    	match bit {
			Some(0) => "A".to_string(), //0b00
			Some(1) => "C".to_string(), // 0b01
			Some(2) => "G".to_string(), // 0b10
			Some(3) => "T".to_string(), // 0b11
			Some(_) => panic!("error in decoding binary nucl!"),
			None => "n".to_string(),
		}
    }


    /// splits the cigar string into CigarTuples. They are a representation of the (\d+)([MIDX]),
    /// but also store their position in both the Vec<CigarEnum> and the cigar.cigar string.
	pub fn to_cigar_tupel_vec( &self, only_gaps:bool ) -> Vec<CigarTuple> {
		self.str_to_tuple_vec( &self.cigar, self.state_changes, only_gaps)
	}

    fn str_to_tuple_vec( &self, cig:&str, state_changes:usize, only_gaps:bool ) -> Vec<CigarTuple> {

		let re = CigarEnum::get_regex(); // Example CIGAR regex
		let mut cigar_tuples = Vec::with_capacity( state_changes );
		let mut vec_pos = 0;
		let mut str_pos = 0;
		let mut read_pos = 0;
		let mut database_pos = 0;

	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter( cig ) {

			let length: usize = cap[1].parse().unwrap();
            let cigar_tuple = CigarTuple::from_match( &cap, vec_pos, str_pos, read_pos, database_pos );

            //cigar_tuple.print_debug();
            vec_pos += cigar_tuple.vec_len;
            str_pos += cigar_tuple.str_len;
            if cigar_tuple.option == CigarEnum::Insertion{
            	read_pos += length;
            }else if cigar_tuple.option == CigarEnum::Deletion {
            	database_pos += length;
            }else{
            	read_pos += length;
            	database_pos += length;
            }
            // store the value if we need it
            if  only_gaps && cigar_tuple.option.is_gap()  {
                cigar_tuples.push(cigar_tuple);
            }else if ! only_gaps{
            	cigar_tuples.push(cigar_tuple);
            }
	    }
	    cigar_tuples
	}


	pub fn to_vec (&self) ->  Vec<CigarEnum> {
		return self.string_to_vec( &self.cigar )
	}


    pub fn string_to_vec(&self, cigar_string:&str ) -> Vec<CigarEnum> {
    	let mut result = Vec::with_capacity( self.len() );
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    	for cap in re.captures_iter(cigar_string) {
	        let count: usize = cap[1].parse().unwrap();
	        let operation = &cap[2];
	        let cigar_type = match CigarEnum::from_str(operation) {
	        	None => { 
	        		panic!("You try to revert the Cigar String {cigar_string} back to a vector of CigarEnum - but we do not recognize the {operation} as a CigarEnum!")
	        	},
	            Some (cig) => {
	            	cig
	            },
	        };

	        // Add the operation to the result vector `count` times
	        result.extend(std::iter::repeat(cigar_type).take(count));

   		}
   		result
    }

    pub fn restart_from_cigar(&mut self, cigar_string:&str ) {
    	self.clear();
    	self.cigar = cigar_string.to_string();

    	let result = self.to_cigar_tupel_vec( false );

	    self.reset_fom_path( &result );
    }

    	/// This function returns the number of equal bases before the position and therefore
	/// needs both sequences and the positions to start from on both elements.
	fn neg_look_ahead<T>( &self, read:&T, database:&T, pos_r: usize, pos_d:usize ) -> usize
	where
    T: BinaryMatcher{
		let mut pos_read = pos_r;
		let mut pos_db = pos_d;
		let mut overlaps = 0;
		let mut match_a = "".to_string();
		let mut match_b = "".to_string();
 		while read.get_nucleotide_2bit( pos_read ) ==  database.get_nucleotide_2bit(  pos_db ) {
			overlaps +=1;
			match_a += &self.get_nucleotide_2bit(read.get_nucleotide_2bit( pos_read )  );
			match_b += &self.get_nucleotide_2bit(database.get_nucleotide_2bit( pos_db )  );
			if pos_read == 0 || pos_db == 0 {
				break;
			}
			pos_read -= 1;
			pos_db -=1;
		}
		#[cfg(all(debug_assertions, feature = "detailed_mapping_debug"))]
		{
			println!("searching for neg_look_ahead from pos_r {pos_r} and pos_d {pos_d} found {overlaps} ({match_a} and {match_b}) overlaps ({pos_read} and {pos_db} were my iterators");
			println!("that is the mismatching pair: {}, {}", &self.get_nucleotide_2bit(read.get_nucleotide_2bit( pos_read )),  &self.get_nucleotide_2bit(database.get_nucleotide_2bit( pos_db ) ));
		}
		overlaps
	}




    /// cleanup in this regards is meant to fix alignment errors.
    pub fn clean_up_cigar<T>(&mut self, seq1:&T, seq2:&T)
    where
    	T:BinaryMatcher + std::fmt::Display,
    {
    	if self.fixed.is_some() {
    		//eprintln!("This cigar was already fixed!");
   			return ();
    	}
    	//( self.dropped_start, self.dropped_end ) = seq1.get_dropped_values();
    	self.fix_di_problems( 0, seq1, seq2 );
        //println!("Before the soft clip I have this result {self} for these sequences:{seq1}\n {seq2}\n");
        self.soft_clip_start_end();

	    //println!("After the soft clip I have this result {self} for these sequences:{seq1}\n {seq2}\n");
    }

    pub fn mapping_quality(&self ) -> u8 {

    	if &self.cigar == ""{
    		return 0
    	}
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    
	    // Initialize counts for 'M' and other operations
	    let mut m_count = 0;
	    let mut other_count = 0;
	    
	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter(&self.cigar) {
	        let count: usize = cap[1].parse().unwrap();
	        let operation = &cap[2];
	        match operation {
	            "M" => m_count += count,
	            //"S" => (),// soft clipped is ignored here
	            //"H" => (),// hard clipped is ignored here
	            "N" => (),// not matched (intron) is ignored, too
	            _ => other_count += count,
	        }
	    }
	    if other_count == 0 {
	    	//println!("MapQ = {}",40);
	    	return 40_u8
	    }
	    let ratio = (m_count as f32) / ((other_count + m_count) as f32 );
	    //println!("MapQ = {}",  (40.0  * ratio )as u8 );
	    (40.0  * ratio ) as u8
    }

    pub fn edit_distance(&self ) -> f32 {
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    
	    // Initialize counts for 'M' and other operations
	    let mut other_count = 0.0;
	    let mut total = 0.0;	    
	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter(&self.cigar) {
	        let count: f32 = cap[1].parse().unwrap();
	        total += count;
	        let operation = &cap[2];
	        match operation {
	            "M" => (), // match
	            //"S" => total -= count,// soft clipped is igniored here
	            //"H" => total -= count,// hard clipped is igniored here
	            "N" => total -= count,// not matched (intron) also ignored
	            "=" => (), // not quite sure, but should likely be good too - or?
	            _ => other_count += count,
	        }
	    }
	    other_count  / total
    }

    pub fn mapped(&self) -> f32{
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    	let mut total = 0.0;
    	for cap in re.captures_iter(&self.cigar) {
    		let count: f32 = cap[1].parse().unwrap();
    		match &cap[2]{
    			"M" => total += count,
    			_ => (),//only count MATCH
    		}
    	}
    	total
    }

    pub fn as_alignement<T: BinaryMatcher>(&self, read:&T, database:&T) ->String {
    	let mut cigar = "".to_string();
    	let mut tens =  "".to_string();
    	let mut ones = "".to_string();
    	let mut a = "".to_string();
    	let mut b = "".to_string();
    	let mut pos_a = 0;
    	let mut pos_b = 0;
    	let mut id = 0;

    	for operation in self.string_to_vec( &self.cigar ) {
    		cigar += &operation.to_string();
    		tens += &(id / 10).to_string();
    		ones += &(id % 10).to_string();
    		id += 1;
    		if id ==100 {
    			id = 0;
    		}
    		match operation{
    			CigarEnum::Match => {
    				let nucl_a = self.get_nucleotide_2bit( read.get_nucleotide_2bit(pos_a) );
					let nucl_b = self.get_nucleotide_2bit( database.get_nucleotide_2bit(pos_b) );
    				a += &nucl_a;
    				b += &nucl_b;
    				pos_a +=1;
    				pos_b +=1;
    			},
				CigarEnum::Mismatch => {
					//"X",
					let nucl_a = self.get_nucleotide_2bit( read.get_nucleotide_2bit(pos_a) );
					let nucl_b = self.get_nucleotide_2bit( database.get_nucleotide_2bit(pos_b) );
    				a += &nucl_a;
    				b += &nucl_b;
    				pos_a +=1;
    				pos_b +=1;
				},
				CigarEnum::Insertion | CigarEnum::Hardclip | CigarEnum::Softclip => {
					//"I",
					let nucl_a = self.get_nucleotide_2bit( read.get_nucleotide_2bit(pos_a) );
					a += &nucl_a;
					b += "-";
					pos_a +=1;
				},
				CigarEnum::Deletion | CigarEnum::Nothing => {
					//"D",
					let nucl_b = self.get_nucleotide_2bit( database.get_nucleotide_2bit(pos_b) );
					b += &nucl_b;
					pos_b +=1;
					a +="-";
				},
				CigarEnum::Empty => panic!("There is an empty cigar entry in your vector!"),
    		}
    	}
    	format!("Cigar String: {}\n{}\n{}\n{}\n{}\n{}",self.cigar, tens, ones, a, b , cigar )
    }

    pub fn check_alignment<T: BinaryMatcher>(&mut self, read:&T, database:&T) {

    	let mut changed = false;
    	let mut enum_vec = self.string_to_vec( &self.cigar );
    	let mut pos_a = 0;
    	let mut pos_b = 0;
    	for operation in enum_vec.iter_mut()  {

    		match operation{
    			CigarEnum::Match => {
					if self.get_nucleotide_2bit( read.get_nucleotide_2bit(pos_a) ) 
						!= 
						self.get_nucleotide_2bit( database.get_nucleotide_2bit(pos_b) ) {
						*operation = CigarEnum::Mismatch;
						changed = true;
					}
					pos_a +=1;
					pos_b +=1;
    			},
				CigarEnum::Mismatch => {
					//"X",
					if self.get_nucleotide_2bit( read.get_nucleotide_2bit(pos_a) )
					   ==
					   self.get_nucleotide_2bit( database.get_nucleotide_2bit(pos_b) ) {
						*operation = CigarEnum::Match;
						changed = true;
					}
    				pos_a +=1;
    				pos_b +=1;
				},
				CigarEnum::Insertion | CigarEnum::Hardclip | CigarEnum::Softclip => {
					pos_a += 1; // Insertion advances only read position
				},
				CigarEnum::Deletion | CigarEnum::Nothing => {
					pos_b += 1; // Deletion advances only database position
				},
				CigarEnum::Empty => panic!("There is an empty cigar entry in your vector!"),
    		}
    	}

    	if changed {
    		self.reset_fom_path( &Self::collapse_cigar_enum_vec(&enum_vec) ) ;
    	}
    }

    pub fn finalize(&mut self ) {
    	let mut path = Self::collapse_cigar_enum_vec ( &self.to_vec());
    	// let's see if there are some useless D/I or I/D combinations in the same length
    	// replace them by one nX instead.
    	for id in 0..path.len()-1 {
    		if path[id].is_gap() && path[id].option.opposite( &path[id+1].option) && path[id].len() == path[id +1 ].len() {
    			path[id].option = CigarEnum::Mismatch;
    			path[id+1].vec_len = 0;
    		}
    	}
    	self.reset_fom_path( &path );
    }

    #[allow(dead_code)]
    /// expand a CigarTupe vector into a CigarEnum Vector
    fn expand_cigar_tuple_vec( vec: &[CigarTuple]) -> Vec<CigarEnum> {
    	let mut ret = Vec::<CigarEnum>::with_capacity( vec.into_iter().map(|v| v.len() ). sum::<usize>() );
    	for option in vec {
    		option.extend_vec( &mut ret);
    	}
    	ret
    }

    /// collapse a CigarTupe vector into a CigarEnum Vector
    fn collapse_cigar_enum_vec( vec: &[CigarEnum] ) -> Vec<CigarTuple> {
    	let mut path = Vec::<CigarTuple>::with_capacity( 30 );
    	for option in vec {
    		if let Some(tuple) = path.last_mut() {
			    if tuple.is_a(&option) {
			        tuple.vec_len += 1;
			    } else {
			        path.push(CigarTuple::from_scratch(*option, 1));
			    }
			} else {
				// we obviousely are at position 0:
			    path.push(CigarTuple::from_scratch(*option, 1));
			}
    	}
    	path
    }

    pub fn to_sam_string(&mut self, length:usize ) -> Option<(String,  usize )>
    {

    	self.fixed=None;
    	self.soft_clip_start_end();
    	let mut ret = self.cigar.to_string();

		/*let re_start = Regex::new(r"^(\d+)([ID])").unwrap();

		let move_start = if let Some(mat) =re_start.captures(&ret) {
			let clippable = mat.get(1).unwrap();
			let length : usize= mat.get(1).unwrap().as_str().parse().unwrap();
			let option = mat.get(2).unwrap().as_str();
			match option{
				"I" => {
					// crap - this will be a negative change for the read location!
					ret.replace_range(clippable.start()..clippable.end()+1, &format!("{}S", length) );
					0
				},
				"D" => {
					// move the match n bp instead!
					ret.replace_range(clippable.start()..clippable.end()+1, "" );
					length
				},
				_ => unreachable!()
			}
		}else {
			0
		};

		let re_end = Regex::new(r"(\d+)([DI])$").unwrap();
		if let Some(mat) =re_end.captures(&ret) {
			let clippable = mat.get(1).unwrap();
			let length : i64= mat.get(1).unwrap().as_str().parse().unwrap();
			let option = mat.get(2).unwrap().as_str();

			match option{
				"I" => {
					ret.replace_range(clippable.start()..clippable.end()+1, &format!("{}S", length) );
				},
				"D" => {
					// that can just go.
					ret.replace_range(clippable.start()..clippable.end()+1, "" );
				},
				_ => unreachable!()
			}
		}*/

		// and now we need to handle the possibiliity that the mapper has cut parts of our string, too.
		if self.dropped_start > 0 {
			ret = format!("{}H{}", self.dropped_start, ret);
		}
		if self.dropped_end > 0 {
			ret += &format!("{}H", self.dropped_end);
		}

		// todo("If the start Clips would change the positions?")
		let (mine, _other) = self.calculate_covered_nucleotides( &ret );
		if mine != length {
			// likely a really really crappy mapping anyhow - so just ignore that
		    // eprintln!("The cigar does not have the correct length {}! {}", length, ret );
			return None
		}
		return Some( (ret, 0 ) )
    }

    /// The new mapper likes to add DDJJ elements that are basically
    /// --tt
    /// tt-- combinations. So to say just mapping errors.
    /// They need to go!
    #[allow(unused_variables)]
    pub fn fix_di_problems<T>( &mut self, _mapping_start: usize, read:&T, database:&T )
    where
    T: BinaryMatcher{
    	if self.contains.iter().all( |&t| t ) && self.state_changes > 20 {
    		// too crappy read - ignore
    		return
    	}
    	if self.fixed.is_some() {
    		return;
    	}
    	#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    	{
    		println!("########################\nStarting a new FIX\n########################\n{}",self.as_alignement(read, database));
    		println!("\nFixing fix_di_problems locations in this alignement:\n{}", self.as_alignement(read, database) );
    	}



    	// we look for VERY localized errors - hence we would not create issues prior to the opsition we work on
    	// if we wrok from the end on.
    	if self.contains[CigarEnum::Insertion.to_id()] && self.contains[CigarEnum::Deletion.to_id()]{

    		// we need this Cigar as vec<CigarEnum> as this is easiest to change
    		//let cigar_enum_vec = self.to_vec();

    		// and we need the opsitions of all but not only the gaps:
    		let mut cigar_tuple_vec = self.to_cigar_tupel_vec( false );

    		// would be interesting if we need to regenerate the Cigar string in the end or not.
    		let mut modified = false;
    		if cigar_tuple_vec.len() < 2{
    			return
    		}
    		for id in (1..cigar_tuple_vec.len()).rev() {
    			if
    			// there are never two entries of the same type after each other
    			// and as there are only Insertion and Deletion elements they need to be either one of them
    			cigar_tuple_vec[id].is_gap()
    			// If they are both that means we have a DI or ID combo!
    			&& cigar_tuple_vec[id.saturating_sub(1)].is_gap() 
    			// and the length of both of them needs to be the same
    			&& cigar_tuple_vec[id].len() == cigar_tuple_vec[id.saturating_sub(1)].len() 
    			{
    				//let len = cigar_tuple_vec[id].len();
    				let (on_read, on_db) = match cigar_tuple_vec[id].option{
    					CigarEnum::Insertion => {
    						( 
    							cigar_tuple_vec[id].slice_from_read(read),
    							cigar_tuple_vec[id.saturating_sub(1)].slice_from_database(database),
    							)
    					},
    					CigarEnum::Deletion =>{
    						( 
    							cigar_tuple_vec[id.saturating_sub(1)].slice_from_read(read),
    							cigar_tuple_vec[id].slice_from_database(database),
    							)
    						
    					},
    					_ => {
    						panic!("{} is not a gap!", cigar_tuple_vec[id])
    					}
    				};
    				
    				
    				if on_read == on_db  {
    					//println!("And they even had the same nucleotides");
    					let mut changed = cigar_tuple_vec.remove(id);
    					changed.option = CigarEnum::Match;
    					cigar_tuple_vec[id-1] = changed;
    					modified = true;
    					#[cfg(all(debug_assertions, feature = "mapping_debug"))]
	    				{
	    					let mut tmp = self.clone();
	    					tmp.reset_fom_path( &cigar_tuple_vec );
	    					println!("Foud a DI or ID problem! at read pos {:?} database_pos {:?} - the result:\n{}", on_read, on_db, tmp.as_alignement(read, database) );
	    					cigar_tuple_vec[id].print_debug();
	    				}
    				}else {
    					//println!("the two slices did not match: {:?} vs {:?} :-(",on_read  ,  on_db )
    				}
    				//println!("next one\n##########################\n");
    			}
    		}

    		if modified {
        		let cigar = cigar_tuple_vec.iter()
        			.map(|cigar| cigar.to_string()) // Convert each CigarTuple to a string
        			.collect::<String>(); // Collect into a single string
        		self.restart_from_cigar(&cigar );

        		#[cfg(all(debug_assertions, feature = "mapping_debug"))]
        		println!("The modified alignement: \n{}\n", self.as_alignement(read, database));
    		}
    	}

	    // next step - D an I's tend to be on the wrong side of a gap - creating a ton of extra changes!
	    //------------------------------------------XXXXX-------------------------------------------
	    //GGTGTGACCATGTTCATTATAATCTCAAAGGAGAAAAAAAAACCTT-GTAAAAAAAAGCAAAAACAACAACAAAAAAACAATCTTATTCC
		//GGTGTGACCATGTTCATTATAATCTCAAAGGAGAAAAAAAAAACCTTGTAAAAAAAAGCAAAAACAACAACAAAAAAACAATCTTATTC-
		//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXMXMDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMI
		

		// that could fix some of these issues already...
    	self.fix_next_gap_location( read, database, 0 );
		
		#[cfg(all(debug_assertions, feature = "mapping_debug"))]
		println!("Fixed DI problems locations in this alignement:\n{}\n{self}\n######################################", self.as_alignement(read, database) );
		self.fixed = Some( CigarEndFix::Na );
	}


	/// checks if a deletion or insertion would be beneficial for a (modified) mapping
	/// by checking if the (at least) 2 bp flanking on each side would be mapping with the deletion
	pub fn is_good_gap<T>(&self, gap:&CigarTuple, read:&T, database:&T )->bool 
	where
	T: BinaryMatcher{
		//println!("is_good_gap? read_pos {}; database_pos {}\n{}", gap.read_position, gap.database_position, self.as_alignement( read, database) );

		if gap.option == CigarEnum::Insertion{
			let match_to_start = self.neg_look_ahead(read, database, 
				gap.read_position-1, gap.database_position-1 );
			let match_behind = self.neg_look_ahead(read, database, 
				gap.read_position + gap.len() +1, gap.database_position + 1);

			match_to_start > 1 && match_behind == 2

		}else if gap.option == CigarEnum::Deletion {
			let match_to_start = self.neg_look_ahead(read, database, 
				gap.read_position-1, gap.database_position-1 );
			let match_behind = self.neg_look_ahead(read, database, 
				gap.read_position +1, gap.database_position + gap.len() +1 );

			match_to_start > 1 && match_behind == 2
		}else {
			false
		}
	}

	/// fills in matches entries in cigar_tuple_vec with current_tuple while ignoring not_touch and removing not_touch.opposite()
	/// It starts at id and tranverses the vector in reverse.
	/// returns the amount of entries that would need to be flipped back - if applicable and the id to start at if
	/// reverting is of interest.
	#[allow(unused_variables)]
	fn replace_n_and_start_at<T>(&mut self,  cigar_tuple_vec: &mut Vec<CigarTuple>, current_tuple: CigarEnum, 
		kill: Option<CigarEnum>, to_flip: usize, id:usize, read:&T, database:&T) -> Result<usize,String> where
    T: BinaryMatcher{

		let mut this = id;
		#[allow(unused_variables)]
		let mut dropped = 0;
		let mut added = 0;
		let mut matches = to_flip;
		//let mut last_gap;
		let mut skip = 0;

		let mut flip_back:i32 = 0;
		#[cfg(all(debug_assertions, feature = "mapping_debug"))]
		println!("replace_n_and_start_at introducing {} while removing {:?}",current_tuple, kill);

		while matches > 0{
			//last_gap = matches;
			//println!("Still some way to go: {matches} with this == {this} and flip_back {flip_back}");
			if kill
			    .as_ref()
			    .map_or(false, |nt| nt == &cigar_tuple_vec[this].option)
			{
				//println!("Removoing an opponent: {}",cigar_tuple_vec[this] );

				dropped += cigar_tuple_vec[this].len();
				flip_back -=  (cigar_tuple_vec[this].len()) as i32;
				cigar_tuple_vec.remove( this );

			}
			else if cigar_tuple_vec[this].option.is_gap(){
				//println!("Found the gap {}", cigar_tuple_vec[this] );
				match cigar_tuple_vec[this].split_at( matches, current_tuple){
					Some(new_fragment) => {
						matches -= new_fragment.len();
						flip_back += (new_fragment.len()) as i32;
						cigar_tuple_vec.insert(this+1, new_fragment);
					},
					None => {
						flip_back += (cigar_tuple_vec[this].len()) as i32;
						matches -= cigar_tuple_vec[this].len();
						this -=1;
						skip +=1;
					}
				}
			}else{
				//println!("Found {}", cigar_tuple_vec[this]);
				match cigar_tuple_vec[this].split_at( matches, current_tuple){
					Some(new_fragment) => {
						matches -= new_fragment.len();
						cigar_tuple_vec.insert(this+1, new_fragment);
					},
					None => {
						matches -= cigar_tuple_vec[this].len();
						this -=1;
						skip +=1;
					}
				}
			}
			#[cfg(all(debug_assertions, feature = "detailed_mapping_debug"))]
			{
				let mut cig = "".to_string();
	    		for  entry in &mut *cigar_tuple_vec{
	    			cig += &entry.to_string();
	    		}
	    		let mut report = self.clone();
	    		report.restart_from_cigar( &cig );
	    		println!("             updated: {matches} with this == {this} and flip_back {flip_back}");
				println!("Intermediate alignement at {this}; {}:\n{}", cigar_tuple_vec[this], report.as_alignement( read, database) );
			}
		}
		// so let's see if we need to re-introduce some
		let flip_to =  match kill {
			Some(v) => v.get_opposite(),
			None => return Ok( skip + added ),
		};

		//println!("####################################\nadd the gaps back in\n####################################\n");
		// this creates problems if the next position in the data would be actually a D and weare moving I's like this:
		//Cigar String: 1X1M2X1I1X1M1X1D3I44M3M31M
		//000000000011111111112222222222333333333344444444445555555555666666666677777777778888888888
		//012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
		//ATTATTGT-GGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
		//TTGT-GGGT---GGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
		//XMXXIXMXDIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

		//replace_n_and_start_at introducing M while removing Some(Deletion)
        //     updated: 0 with this == 8 and flip_back 3
		//Intermediate alignement at 8; :
		//Cigar String: 1X1M2X1I1X1M1X1D3M44M3M31M
		//000000000011111111112222222222333333333344444444445555555555666666666677777777778888888888
		//012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
		//ATTATTGT-GGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
		//TTGT-GGGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTGnnn
		//XMXXIXMXDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
		//A not existing intermediate to explain the problem
		//000000000011111111112222222222333333333344444444445555555555666666666677777777778888888888
		//012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
		//ATTATTGT-GGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
		//TTGT-G---GGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
		//XMXXIXMXDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

		while flip_back > 0 {
			//println!("Found {} - trying to flip {flip_back} entries to {flip_to}", cigar_tuple_vec[this]);
			if cigar_tuple_vec[this].option == flip_to {
				//println!("alread a {} - no change {flip_back} (this={this}", flip_to);
				this = this.checked_sub(1).expect(&format!("#1 Sorry I failed to re-introduce the gap!! {}", self.as_alignement( read, database )));
				skip +=1;
			}else if cigar_tuple_vec[this].option.opposite(&flip_to) {
				// here we need to be extremely careful!
				// would that actually benefit the match?!
				if cigar_tuple_vec[this].len() < flip_back as usize && self.is_good_gap( &cigar_tuple_vec[this], read, database ) 
				{
					// a new issue (is solved by the new if case!
					//000000000011111111112222222222333333333344444444445555555555666666666677777777778
					//012345678901234567890123456789012345678901234567890123456789012345678901234567890
					//AGTGCGGCGAAGGGG------TAAAATATAGCTACAAGGTGTCAGAGATGAAGGTGCGGCCCGCCTAGACGGGCCAGGACC
					//GT-GCGGCGAAGGGGTACAAA--TATAGCTACAAGGTGTCAGAGATGAAGGTGCGGCCCGCCTAGACGGGCCAGGACCnnn
					//XXIMMMMMMMMMMMMDDDDDDIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
					//possible intermediate:
					//000000000011111111112222222222333333333344444444445555555555666666666677777777778
					//012345678901234567890123456789012345678901234567890123456789012345678901234567890
					//AGTGCGGCGAAGGGG------TAAAATATAGCTACAAGGTGTCAGAGATGAAGGTGCGGCCCGCCTAGACGGGCCAGGACC
					//GT-GCGGCGAAGGGGTACAAA-----TATAGCTACAAGGTGTCAGAGATGAAGGTGCGGCCCGCCTAGACGGGCCAGGACC
					//XXIMMMMMMMMMMMMDDDDDDIIIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
					// Actually would need to be resolved like this:
					//000000000011111111112222222222333333333344444444445555555555666666666677777777778
					//012345678901234567890123456789012345678901234567890123456789012345678901234567890
					//AGTGCGGCGAAGGGGTA-AAATATAGCTACAAGGTGTCAGAGATGAAGGTGCGGCCCGCCTAGACGGGCCAGGACC
					//GT-GCGGCGAAGGGGTACAAATATAGCTACAAGGTGTCAGAGATGAAGGTGCGGCCCGCCTAGACGGGCCAGGACC
					//XXIMMMMMMMMMMMMMMDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM


					// this deletion would actually benefit the alignement.
					#[cfg(all(debug_assertions, feature = "detailed_mapping_debug"))]
					println!("This {} would actually benefit the alignement! - I need to still flip {} entries", cigar_tuple_vec[this], flip_back,);
					this = this.checked_sub(1).expect(&format!("#2 Sorry I failed to re-introduce the gap!! {}", self.as_alignement( read, database )));
					skip +=1;
					//000000000011111111112222222222333333333344444444445555555555666666666677777777778888888888
					//012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
					//ATTATTGT-GGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
					//TTGT----GGGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG
					//XMXXIIIIDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
				}else {
					// Nay - crap info.
					#[cfg(all(debug_assertions, feature = "detailed_mapping_debug"))]
					println!("This {} can go! {} and {}", cigar_tuple_vec[this], cigar_tuple_vec[this].read_position -1, 
						cigar_tuple_vec[this].database_position  );
					if flip_back > cigar_tuple_vec[this].vec_len as i32{
						flip_back -= cigar_tuple_vec[this].vec_len as i32;
						cigar_tuple_vec[this].vec_len = 0;
					}else {
						// crap - we can not kill all of that, but need to spliot this, too!
						match cigar_tuple_vec[this].split_at( flip_back.try_into().unwrap(), CigarEnum::Empty ){
							Some(new_fragment) => {
								// this is the empty one - so we can savely ignore that one ;-)
								flip_back = 0;
							},
							None => {
								panic!("We tried to remove part of an opposite gap as the last step of our process!\nAnd failed!\n{}", 
									self.as_alignement( read, database)
									);
							}
						}
					}
				}
				
				this = this.checked_sub(1).expect(&format!("3# Sorry I failed to re-introduce the gap!! {}", self.as_alignement( read, database )));
				skip +=1;

				//println!("This is the opposite - right - killing it end reducing flip_back to {flip_back} (this={this}");
			}
			else {
				match cigar_tuple_vec[this].split_at( flip_back.try_into().unwrap() , flip_to){
					Some(new_fragment) => {
						flip_back -= new_fragment.len() as i32;
						cigar_tuple_vec.insert(this+1, new_fragment);
						added +=1;
						//println!("overly large match!");
						if flip_back > 0 {
							this = this.checked_sub(1).expect(&format!("Sorry I failed to re-introduce the gap!! {}", self.as_alignement( read, database )));
							skip +=1;
						}
					},
					None => {
						//println!("Too small match");
						flip_back -= cigar_tuple_vec[this].len() as i32;
						if flip_back > 0 {
							this = this.checked_sub(1).expect(&format!("4# Sorry I failed to re-introduce the gap!! {}", self.as_alignement( read, database )));
							skip +=1;
						}
					}
				}
			}
			
			/*let mut cig = "".to_string();
    		for  entry in &mut *cigar_tuple_vec{
    			cig += &entry.to_string();
    		}
    		let mut report = self.clone();
    		report.restart_from_cigar( &cig );
    		println!("   re-adding updated: {matches} with this == {this} and flip_back {flip_back}");
			println!("Intermediate alignement at {this}; {}:\n{}", cigar_tuple_vec[this], report.as_alignement( read, database) );
			*/
			if flip_back == 0 {
				break;
			}
		}
		
		self.reset_fom_path( &cigar_tuple_vec );
		//self.check_alignment( read, database);

		#[cfg(all(debug_assertions, feature = "detailed_mapping_debug"))]
		println!("####################################\nAfter the fix I return {} and got the alignement:\n{}\n####################################",skip + added +1, self.as_alignement(read, database));
		
		return Ok( skip + added +1 );
	}

	/// we need to update the read and database positions for this vector
	fn update_positions( &self, vec: &mut Vec<CigarTuple> ) {
		let mut pos_read = 0;
		let mut pos_database = 0;

		for item in vec {
			item.read_position = pos_read;
			item.database_position = pos_database;
			if item.option == CigarEnum::Insertion {
				pos_read += item.len();
			}else if item.option == CigarEnum::Deletion {
				pos_database += item.len();
			}else {
				pos_read += item.len();
				pos_database += item.len();
			}
		}
	}

	/// a streamlined function that processes the match from the rear and checks both deletions and insertions.
	/// Tuns multiple times untill there are no more entries to be checked or there are no more deletions/insertions in the data.
	fn fix_next_gap_location<T>( &mut self, read:&T, database:&T, skip: usize )
	where
    T: BinaryMatcher{
    	if self.state_changes > 30{
    		//useless crap!
    		return
    	}

    	if self.contains[CigarEnum::Deletion.to_id()] || self.contains[CigarEnum::Insertion.to_id()]{
    		//println!("Fixing {} locations in this alignement:\n{}", this_option, self.as_alignement(read, database) );
    		//let mut cigar_vec = self.to_vec();
    		//let mut skipped = skip;
    		let mut fixes = 0;
    		#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    		println!("Starting alignement:\n{}", self.as_alignement( read, database) );
    		let mut cigar_tuple_vec = self.to_cigar_tupel_vec( false );
    		if skip >= cigar_tuple_vec.len() -2 {
    			return;
    		}
    		for id in (1..cigar_tuple_vec.len()-1).rev().skip( skip ) {
				//skipped +=1;
    			if ! cigar_tuple_vec[id].is_gap(){
    				continue;
    			}else if cigar_tuple_vec[id].len() == 0{
    				continue;
    			}else if cigar_tuple_vec[id].database_position == 0 {
    				// this should not be worked on!
    				// need to change the lib!
    				continue;
    			}
    			// make sure all the positions are still up to date!
    			if fixes > 0 {
    				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    				println!("I think it is more save to restart this here! \n{}", self.as_alignement(read, database));
    				// clean out the likely messy inserts / deletions
    				self.finalize();
    				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    				println!("after the finalze run:\n{}", self.as_alignement(read, database));
    				return self.fix_next_gap_location( read, database, 0 );
    			}else {
    				self.update_positions( &mut cigar_tuple_vec );
    			}
    			
    			let current_tuple = cigar_tuple_vec[id].clone();

    			let (on_read, on_db ) = match current_tuple.option{
    				CigarEnum::Insertion => {
    					( current_tuple.read_position + current_tuple.len() - 1 , current_tuple.database_position - 1 )
    				},
    				CigarEnum::Deletion =>  {
    					( current_tuple.read_position -1, current_tuple.database_position + current_tuple.len() -1 )
    				},
    				_ => {
    					panic!("This is no gap enum: {}", current_tuple.option );
    				}

    			};

    			let matches = self.neg_look_ahead(read, database, on_read, on_db );
    			if matches == 0 {
    				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    				println!("no overlap for test {on_read} {on_db}");
    				continue;
    			}
    			if cigar_tuple_vec[id-1].option.opposite( &current_tuple.option ) 
    				&& matches <= current_tuple.len()
    				&& matches <= cigar_tuple_vec[id-1].len()
    				{
    				// this is by far not worth looking into!
    				cigar_tuple_vec[id].vec_len -= matches;
    				cigar_tuple_vec[id-1].vec_len -= matches;
    				continue;
    			}
				fixes +=1;
    			#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    			{
    				println!("id:{}\n{}\nI found a gap {} at read {} db {} and the overlap was {}",
    				 id,
	    			 self.as_alignement( read, database),
	    			 current_tuple, 
	    			 on_read, 
	    			 on_db,
	    			 matches );
	    			current_tuple.print_debug();
	    		}

	    		// probably better to keep that on the level of touples again.

	    		let _skip_more = match self.replace_n_and_start_at( &mut cigar_tuple_vec, CigarEnum::Match, 
	    			Some(current_tuple.option.get_opposite()), matches, id, read, database ){
	    			Ok(ret) => ret,
	    			Err(e) => {
	    				panic!("replace_n_and_start_at to {} evade {:?}: {e}\n{} and flip {matches}",  CigarEnum::Match, Some(current_tuple.option), self.as_alignement( read, database) );
	    			}
	    		};

    			#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    			println!("After the replace_n_and_start_at run we will go on on id {}-1:\n{}\n####################################\n",
    				id, self.as_alignement( read, database));
    			if ! self.contains[CigarEnum::Mismatch.to_id()]{
    				break;
    			}
    		}

		}
	}


    /// This will soft clip (or better replace with X's) unstable alignements that shift between states frequently and have a low matching count
    /// These areas are identified by those RegExp constructs:
    /// let start = r"^((?:[1-7]M|[1-9][0-9]*[IXD]){4,})";
    /// let end = r"((?:[1-7]M|[1-9][0-9]*[IXD]){4,})$";
    pub fn soft_clip_start_end( &mut self) {
    	// actually that is a REALLY bad idea!

    	if self.state_changes < 6 && self.fixed == None{
    		self.fixed = Some(CigarEndFix::Na);
    		return;
    	}
    	if self.fixed != None{
    		//println!("Do not run soft_clip_start_end twice! {:?}", self.fixed );
    		return;
    	}
    	let start = r"^((?:[1-7]M|[1-9][0-9]*[IXD]){4,})";
    	//let start = r"^((?:[123456789][IXMD]){5,})";
    	let re_start = regex::Regex::new(start).unwrap();

		let end = r"((?:[1-7]M|[1-9][0-9]*[IXD]){4,})$";
    	//let end = r"((?:[123456789][IXMD]){5,})$";
    	let re_end = regex::Regex::new(end).unwrap();

    	let old_cigar = &self.cigar.clone();
		self.fixed = Some(CigarEndFix::Na);

		//println!("Just while debugging: the 'old cigar' = {old_cigar} and self = {self}");

	    if let Some(mat) = re_end.captures(&old_cigar) {
	        if let Some(clippable) = mat.get(1) {
	        	if clippable.start() > 0{
	        		let prev_char:&char = &self.cigar[(clippable.start()-1)..clippable.start()].chars().next().unwrap();
				    if prev_char.is_digit(10) {
				        // The character is a digit - overmatched
				        let clipped_part = &self.cigar[(clippable.start()+2)..clippable.end()];

						#[cfg(debug_assertions)]	        
				        let (mine, _other) = self.calculate_covered_nucleotides(clipped_part);
				        #[cfg(not(debug_assertions))]
				        let (mine, _other) = self.calculate_covered_nucleotides(clipped_part);

				        self.cigar.replace_range((clippable.start()+2)..clippable.end(), &format!("{}X", mine));
				        self.fixed = Some(CigarEndFix::End);
				    } else {
				        // The character is not a digit - great
				        let clipped_part = &self.cigar[clippable.start()..clippable.end()];
		        		let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
		        		//println!("I'll clip this part: {clipped_part} with a mine length of {mine}");
		            	self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}X", mine));
		            	self.fixed = Some(CigarEndFix::End);
				    }
					
				}else {
			        // The character is not a digit - great
			        let clipped_part = &self.cigar[clippable.start()..clippable.end()];
			        let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
			        //println!("not an digit - I'll clip this part: {clipped_part} with a mine length of {mine}");
			        self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}X", mine));
			        self.fixed = Some(CigarEndFix::End);
			    }
			    
	        }
	    }
		#[cfg(debug_assertions)]
		println!("Updated cigar after end clipping: {self}");

	    let problem = r"^\d\d*S$";
		let re_problem = regex::Regex::new(problem).unwrap();

	   	if let Some(_mat) = re_problem.captures(&self.cigar.clone()) {
	   		//panic!("This should not happen in the tests!");
	   		self.state_changes = 1;
	   		return;
	   	}
	    
	    if let Some(mat) = re_start.captures(&old_cigar) {
	        if let Some(clippable) = mat.get(1) {
	            if clippable.start() <= clippable.end() && clippable.end() <= self.cigar.len() {
	        		let clipped_part = &self.cigar[clippable.start()..clippable.end()];
	        		let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
	            	self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}X", mine));
	            	match self.fixed{
	            		Some(CigarEndFix::End) => self.fixed = Some( CigarEndFix::Both),
	            		Some(CigarEndFix::Na) => self.fixed =Some( CigarEndFix::Start),
	            		_ => unreachable!() ,
	            	}
	            	
	        	}else {
	        		// This can be totally normal for really really crappy matches.
	        		eprintln!("With {self} I found a crappy match {} {}: {} old: {}", clippable.start(), clippable.end(), self.cigar.len(), old_cigar);
	        	}
	            
	        }
	    }
	    #[cfg(debug_assertions)]
	    println!("Updated cigar after start clipping: {self}");

	    self.state_changes = self.cigar.chars().filter(|c| !c.is_digit(10)).count();

	    /*println!("This is the final cigar inside function: {self}");
	    println!("And more specifically the fixed field: {:?}", self.fixed);*/

    }

    /// Soft clip the sequences at start or end of a read that map with a horrible mapping like 1X1M2X1M4X3M1X1D1X1M3X1M2X1M2X1M1X1I3X59M,
    /// instead of that string I would like to see a 23S59M 

	/// calculates the nucleotides on both mine and the other sequence that has passed at the end of the cigar string
	pub fn calculate_covered_nucleotides(&self, cigar_string: &str) -> (usize, usize) {
	    let mut mine = 0;
	    let mut other = 0;
	    let mut current_number = String::new();
	    let mut inserts = 0;
	    let mut deletions = 0;
	    
	    for c in cigar_string.chars() {
	        if c.is_digit(10) {
	            // If the character is a digit, append it to the current number
	            current_number.push(c);
	        } else {
	            // If the character is not a digit, process the operation
	            let count = current_number.parse::<usize>().unwrap_or(1); // Parse the count, default to 1 if parsing fails

	            match c {
	                'M' | '=' | 'X' => {
	                	mine  += count; // Match, mismatch, or sequence match
	                	other += count; 
	                },
	                'I' => {
	                	inserts += count;
	                	//mine += count; // Insertion
	                },
	                'D' => {
	                	deletions += count;
	                	//other += count;// Deletion or intron
	                },
	                'S' => {
	                	mine += count;// Deletion or intron
	                },
	                'H' => {
	                	mine += count;// hard klipped
	                },
	                _ => {}, // Other CIGAR operations (e.g., P)
	            }
	            current_number.clear(); // Clear the current number for the next operation
	        }
	    }   
	    mine = mine + inserts;//.saturating_sub(deletions);
	    other = other + deletions;//.saturating_sub(inserts) ;
	    (mine, other)
	}

	pub fn len (&self) -> usize{
		let mut mine = 0;
		let mut current_number = String::new();
		let mut inserts = 0;
	    //let mut deletions = 0;

		for c in self.cigar.chars() {
	        if c.is_digit(10) {
	            // If the character is a digit, append it to the current number
	            current_number.push(c);
	        } else {
	            // If the character is not a digit, process the operation
	            let count = current_number.parse::<usize>().unwrap_or(1); // Parse the count, default to 1 if parsing fails

	            match c {
	                'M' | '=' | 'X' => {
	                	mine  += count; // Match, mismatch, or sequence match
	                },
	                'I' => {
	                	inserts += count;
	                	//mine += count; // Insertion
	                },
	                'D' => {
	                	//deletions += count;
	                	//other += count;// Deletion or intron
	                },
	                'S' => {
	                	mine += count;// soft klipped
	                },
	                'H' => {
	                	mine += count;// hard klipped
	                }
	                _ => {}, // Other CIGAR operations (e.g., S, H, P)
	            }
	            current_number.clear(); // Clear the current number for the next operation
	        }
	    }   
	    //mine + deletions.saturating_sub(inserts)
	    mine + inserts//.saturating_sub(deletions)
	}



/*
	fn populate_contains( &mut self, cigar: &[CigarEnum] ) {
		for entry in cigar{
			self.contains[entry.to_id()] = true;
		}
	}

	/// this function literally checks for 1D1J or vice versa. An artifact from earlier fixes.
	pub fn fix_1d1i_1i1d(&mut self, cigar: &mut Vec<CigarEnum>, start_pos: Option<usize>) {

		// for a test - just not do this:
		// return;
		if ! self.contains.iter().any(|&x| x){
			self.populate_contains(&cigar);
		}
		if ! self.contains[CigarEnum::Deletion.to_id()] &&  ! self.contains[CigarEnum::Insertion.to_id()] {
			// not necessary to process this!
			return ();
		}

	    let mut start_index = start_pos.clone().unwrap_or(cigar.len() - 1);

	    if start_index >= 3 {
	        for i in 0..=start_index - 3 {
	        	if i > cigar.len() -4 {
	        		break
	        	}
	            if ((cigar[i] == CigarEnum::Match || cigar[i] == CigarEnum::Mismatch) &&
	               (cigar[i + 1] == CigarEnum::Insertion || cigar[i + 1] == CigarEnum::Deletion) &&
	               (cigar[i + 2] == CigarEnum::Match || cigar[i + 2] == CigarEnum::Mismatch) &&
	               (cigar[i + 3] == CigarEnum::Deletion || cigar[i + 3] == CigarEnum::Insertion)) &&
	               cigar[i + 1] != cigar[i + 3] {
	                   
	                   //println!("Found a pattern - {} to {} {}{}{}{}", i + 3, i, cigar[i], cigar[i + 1], cigar[i + 2], cigar[i + 3]);
	                   
	                   // Resolve the pattern
	                   //cigar[i] = CigarEnum::Mismatch;
	                   cigar[i+2] = CigarEnum::Mismatch;
					   cigar[i + 1] = CigarEnum::Mismatch; // each DI ID element represents one X
	                   //println!("Removing {}: {}", i + 3, cigar[i + 3]);
	                   cigar.remove(i + 3);
	                   
	                   //println!("pattern changed to {}{}{}",  cigar[i], cigar[i + 1], cigar[i + 2] );
	                   
	                   // Adjust start_index
	                   start_index += 1;
	            }
	        }
	    }

	    // this sometimes can create 1D1J or 1J1D entries - which in reality are a simple X.
		let mut start_index = start_pos.unwrap_or(cigar.len() - 1);
		if start_index >= 2 {
	        for i in 0..=start_index - 2 {
	        	if i > cigar.len() -2 {
	        		break
	        	}
	        	if 
	               ((cigar[i] == CigarEnum::Insertion || cigar[i] == CigarEnum::Deletion) &&
	               (cigar[i + 1] == CigarEnum::Deletion || cigar[i + 1] == CigarEnum::Insertion)) &&
	               cigar[i] != cigar[i + 1] {
	                   
	                   //println!("Found a pattern - {} to {} {}{}", i , i+1, cigar[i], cigar[i + 1]);
	                   
	                   // Resolve the pattern
	                   cigar[i] = CigarEnum::Mismatch;

	                   //println!("Removing {}: {}", i + 1, cigar[i + 1]);
	                   cigar.remove(i + 1);
	                   
	                   //println!("pattern changed to {}",  cigar[i] );
	                   
	                   // Adjust start_index
	                   start_index += 1;
	            }
	        }
	    }
	}*/



	pub fn calculate_cigar(&mut self, matrix : &Vec<Vec<Cell>>, last_match:bool )  {
			// Trace back the alignment path
		
		    let mut path = Vec::with_capacity( matrix.len().max( matrix[0].len() ));

		    let rows = matrix.len();
		    let cols =  matrix[0].len();

		    let mut i = rows - 1;
		    let mut j = cols - 1;

		    // fix the final entry
		    if last_match {
				path.push(CigarTuple::from_scratch (CigarEnum::Match, 1));
		    }else {
		    	path.push(CigarTuple::from_scratch (CigarEnum::Mismatch, 1));
		    }

		    while i > 0 && j > 0 {
		        let current_score = matrix[i][j].score;
		        //println!("Current score = {current_score}");
		        let diagonal_score = matrix[i - 1][j - 1].score;
		        let up_score = matrix[i - 1][j].score;
		        let left_score = matrix[i][j - 1].score;

		        let max = Self::max3(diagonal_score, up_score, left_score);
		        let option = if max == diagonal_score {
		            i -= 1;
		            j -= 1;
		        	if diagonal_score > current_score {
		        		// mismatch!!!
		        		//println!("adding Mismatch");
		        		CigarEnum::Mismatch
		        	}else {
		        		// match
		        		//println!("adding Match");
		        		CigarEnum::Match
		        	}
		        } else if max == up_score {
		            i -= 1;
		        	//println!("adding Deletion");
		        	CigarEnum::Insertion
		        } else {
		            j -= 1;
		        	//println!("adding Insertion");
		            CigarEnum::Deletion
		        };

		        if let Some(tuple) = path.last_mut() {
				    if tuple.is_a(&option) {
				        tuple.vec_len += 1;
				    } else {
				        path.push(CigarTuple::from_scratch(option, 1));
				    }
				} else {
				    unreachable!();
				}
		    }

		    // If there are remaining gaps at the beginning of the sequences, fill them with the corresponding directions
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
			}

		    if j > 0 {
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

		    path.reverse();

			// Convert the alignment path to CIGAR string
		    
		    //TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCA
		    //TTCATATCGACAATTAGGGTTTACGACCTCGATGTTTCAGGACTAGATAGTA
		    // 37M3D7M3I4M ???! 37M OK 
			self.reset_fom_path(&path);
		    
		    /*
		    let mut csv_table = String::new();

		    for i in 0..rows {
		    	let mut line = "".to_string(); 
		        for j in 0..cols {
		        	line += &format!("{}//{}\t", matrix[i][j].score, matrix[i][j].direction );   
		        }
		        csv_table.push_str(&format!("{}\n", line));
		    }

		    // Write CSV table to a file
		    let mut file = File::create(&format!("{}.csv", &self.cigar)).unwrap();
		    file.write_all(csv_table.as_bytes()).unwrap();

		    println!("You can check my alignement table : '{}'", format!("{}.csv",  &self.cigar ) );
			*/
	}
}