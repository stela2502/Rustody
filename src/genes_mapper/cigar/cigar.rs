//gene_mapper::cigar.rs

use crate::traits::Cell;
use crate::traits::BinaryMatcher;
use crate::genes_mapper::gene_data::GeneData;
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
	dropped_start:usize,
	/// keep track of how many bp have been sliced from this entry's end
	dropped_end:usize,

}

// Implementing Display trait for SecondSeq
impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Cigar {} - fixed {:?}\nstate changed {} - contains {:?}", self.cigar, self.fixed, self.state_changes, self.contains  )
    }
}

impl Default for Cigar {
    fn default() -> Self {
        Cigar {
            cigar: "".to_string(),
            fixed: None,
            debug: false,
            contains: vec![false;4],
            state_changes: 1000,
            dropped_start: 0,
            dropped_end:0,        }
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
		Self{
			cigar: cigar.to_owned(),
			fixed:None,
			debug:false,
			contains: vec![false;4],
			state_changes: 0,
			dropped_start: 0,
            dropped_end:0,
		}
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
		self.contains = vec![false;4];
		self.state_changes = 0;
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
	fn to_cigar_tupel_vec( &self, only_gaps:bool ) -> Vec<CigarTuple> {

		let re = Regex::new(r"(\d+)([MIDX])").unwrap(); // Example CIGAR regex
		let mut cigar_tuples = Vec::with_capacity( self.state_changes );
		let mut vec_pos = 0;
		let mut str_pos = 0;
		let mut read_pos = 0;
		let mut database_pos = 0;

	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter(&self.cigar) {

			let length: usize = cap[1].parse().unwrap();
            let cigar_tuple = CigarTuple::from_match( &cap, vec_pos, str_pos, read_pos, database_pos );

            cigar_tuple.print_debug();
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


	fn to_vec (&self) ->  Vec<CigarEnum> {
		return self.string_to_vec( &self.cigar )
	}


    fn string_to_vec(&self, cigar_string:&str ) -> Vec<CigarEnum> {
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

    	let result = self.string_to_vec(cigar_string);

	    self.reset_fom_path( &result );
	    self.soft_clip_start_end();
    }

    	/// This function returns the number of equal bases before the position and therefore
	/// needs both sequences and the positions to start from on both elements.
	fn neg_look_ahead<T>( &self, read:&T, database:&T, pos_r: usize, pos_d:usize ) -> usize
	where
    T: BinaryMatcher{
		let mut ret = 0;
		while  read.get_nucleotide_2bit( pos_r - ret ) ==  database.get_nucleotide_2bit( pos_d - ret ) {
			ret += 1;
		}
		//println!("searching for neg_look_ahead from pos_r {pos_r} and pos_d {pos_d} found {ret} overlaps");
		ret
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
    	( self.dropped_start, self.dropped_end ) = seq2.get_dropped_values();
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
	    	println!("MapQ = {}",40);
	    	return 40_u8
	    }
	    let ratio = (m_count as f32) / ((other_count + m_count) as f32 );
	    println!("MapQ = {}",  (40.0  * ratio )as u8 );
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
				CigarEnum::Insertion => {
					//"I",
					let nucl_a = self.get_nucleotide_2bit( read.get_nucleotide_2bit(pos_a) );
					a += &nucl_a;
					b += "-";
					pos_a +=1;
				},
				CigarEnum::Deletion => {
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

    pub fn to_sam_string(&self) -> (String,  usize ){
    	let mut ret = self.cigar.to_string();

		let re_start = Regex::new(r"^(\d+)([ID])").unwrap();

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
		}

		// and now we need to handle the possibiliity that the mapper has cut parts of our string, too.
		if self.dropped_start > 0 {
			ret = format!("{}H{}", self.dropped_start, ret);
		}
		if self.dropped_end > 0 {
			ret += &format!("{}H", self.dropped_end);
		}

		return (ret, move_start )
    }

    /// The new mapper likes to add DDJJ elements that are basically
    /// --tt
    /// tt-- combinations. So to say just mapping errors.
    /// They need to go!
    pub fn fix_di_problems<T>( &mut self, mapping_start: usize, read:&T, database:&T )
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
    		println!("Fixing fix_di_problems locations in this alignement:\n{}", self.as_alignement(read, database) );
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
    			&& cigar_tuple_vec[id-1].is_gap() 
    			// and the length of both of them needs to be the same
    			&& cigar_tuple_vec[id].len() == cigar_tuple_vec[id-1].len() 
    			{
    				let len = cigar_tuple_vec[id].len();

    				let (on_read, on_db) = match cigar_tuple_vec[id].option{
    					CigarEnum::Insertion => {
    						( 
    							cigar_tuple_vec[id].slice_from_read(read),
    							cigar_tuple_vec[id-1].slice_from_database(database),
    							)
    					},
    					CigarEnum::Deletion =>{
    						( 
    							cigar_tuple_vec[id-1].slice_from_read(read),
    							cigar_tuple_vec[id].slice_from_database(database),
    							)
    						
    					},
    					_ => {
    						panic!("{} is not a gap!", cigar_tuple_vec[id])
    					}
    				};
    				//println!("Foud a DI or ID problem! at read pos {:?} database_pos {:?}", on_read, on_db );
    				cigar_tuple_vec[id].print_debug();

    				if on_read == on_db  {
    					//println!("And they even had the same nucleotides");
    					let mut changed = cigar_tuple_vec.remove(id);
    					changed.option = CigarEnum::Match;
    					cigar_tuple_vec[id-1] = changed;
    					modified = true;
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
        		//println!("The modified alignement: {}", self.as_alignement(read, database));
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
		println!("Fixed DI problems locations in this alignement:\n{}\n{self}", self.as_alignement(read, database) );
		self.fixed = Some( CigarEndFix::Na );
	}



	/// a streamlined function that processes the match from the rear and checks both deletions and insertions.
	/// Tuns multiple times untill there are no more entries to be checked or there are no more deletions/insertions in the data.
	fn fix_next_gap_location<T>( &mut self, read:&T, database:&T, skip: usize )
	where
    T: BinaryMatcher{

    	if self.contains[CigarEnum::Deletion.to_id()] || self.contains[CigarEnum::Insertion.to_id()]{
    		//println!("Fixing {} locations in this alignement:\n{}", this_option, self.as_alignement(read, database) );
    		let mut cigar_vec = self.to_vec();
    		let mut skipped = skip;

    		let mut cigar_tuple_vec = self.to_cigar_tupel_vec( false );
    		if skip >= cigar_tuple_vec.len() -2 {
    			return;
    		}
    		for id in (1..cigar_tuple_vec.len()-1).rev() {
				skipped +=1;
    			if ! cigar_tuple_vec[id].is_gap(){
    				continue;
    			}
    			let mut current_tuple = cigar_tuple_vec[id].clone();
    			let inv_option = current_tuple.option.get_opposite();

    			let (on_read, on_db ) = match current_tuple.option{
    				CigarEnum::Insertion => {
    					( current_tuple.read_position + current_tuple.len() -1 , current_tuple.database_position-1 )
    				},
    				CigarEnum::Deletion =>  {
    					( current_tuple.read_position-1, current_tuple.database_position + current_tuple.len()-1 )
    				},
    				_ => {
    					panic!("This is no gap enum: {}", current_tuple.option );
    				}

    			};

    			let mut matches = self.neg_look_ahead(read, database, on_read, on_db );

    			//#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    			{
    				println!("{}\nI found a gap {} at read {} db {} and the overlap was {}",
	    			 self.as_alignement( read, database),
	    			 current_tuple, 
	    			 on_read, 
	    			 on_db,
	    			 matches );
	    			current_tuple.print_debug();
	    		}

	    		// probably better to keep that on the level of touples again.

	    		// see if we can finish this element
	    		if cigar_tuple_vec[id].len() < matches {
	    			// just use this as my matches but we have sucked up to much of one sequence
	    			println!("{} can not suck up out matches {matches}", cigar_tuple_vec[id]);
	    			matches -= cigar_tuple_vec[id].len();
	    			// we only need to revert the additinal numcleotide's matches
	    			cigar_tuple_vec[id].option = CigarEnum::Match
	    		}else {
	    			// Oh good - we are finished here!
	    			println!("This {} is more space than we need - finish it up here", cigar_tuple_vec[id]);
	    			let mut temp = cigar_tuple_vec[id].clone();
	    			temp.vec_len -= matches;
	    			cigar_tuple_vec.insert( id, temp);
	    			cigar_tuple_vec[id+1] = CigarTuple::from_scratch( CigarEnum::Match, matches );
	    			println!("We have split the data up into {} and {} ", cigar_tuple_vec[id], cigar_tuple_vec[id+1]);
	    			skipped +=1;
	    			matches = 0;
	    		}
	    		let mut this = id;
	    		println!("Still some way to go: {matches}");
    			while matches > 0 {
    				if this == 0 {
    					panic!("Lib error - I could not move my gap!\n{}", self.as_alignement(read, database));
    				}
    				skipped +=1;
    				if current_tuple.option.opposite( &cigar_tuple_vec[this].option ){
    					println!("Removoing an opponent: {}",cigar_tuple_vec[this] );
    					cigar_tuple_vec.remove( this );
    					continue;
    				}
    				if cigar_tuple_vec[this].option == current_tuple.option {
    					// this can suck up the rest of my info!

    					println!("Found the same gap - but I need to convert others :-(!");
    					this -=1;
    					continue;
    				}
    				if cigar_tuple_vec[this].len() < matches {
    					println!("Found {} - I need MOORE :-(!", cigar_tuple_vec[this]);
    					cigar_tuple_vec[this].option = current_tuple.option;
    					matches -= cigar_tuple_vec[this].len();
    					this -=1;
    				}else {
    					// now we hit a tuple with more entries than we need which is not a gap
    					// create a new tupel and add that behind the one we have at hand
    					println!("Found {} - this is all I had needed...", cigar_tuple_vec[this]);
    					let temp = CigarTuple::from_scratch(current_tuple.option , matches );
    					if this + 1 < cigar_tuple_vec.len() {
				            cigar_tuple_vec.insert(this + 1, temp);
				        } else {
				            cigar_tuple_vec.push(temp); // Push if at the end
				        }
				        skipped +=1;
    					cigar_tuple_vec[this].vec_len -= matches;
    					matches = 0;
    				}
    			}

    			let mut cig = "".to_string();
    			for entry in &cigar_tuple_vec{
    				cig += &entry.to_string();
    			}
    			self.restart_from_cigar( &cig );
    			//#[cfg(all(debug_assertions, feature = "mapping_debug"))]
    			println!("After the fix:\n{}", self.as_alignement( read, database));
 				
    		}

    		println!( "Skip and skipped: {} {}", skip, skipped);
    		if skipped == skip || skip > self.state_changes {
    			return
    			
    		}
    		// there were more changes to be checked
    		self.fix_next_gap_location( read, database, skipped )
		}
	}

    pub fn fix_border_insertion( &mut self, mapping_start: usize, read:&GeneData, database:&GeneData )->usize{
    	return 0;
    	let mut ret =0 ;
    	if self.contains[CigarEnum::Insertion.to_id()] {
			let re_start = Regex::new(r"^(\d+)I").unwrap();
			if let Some(mat) =re_start.captures(&self.cigar) {
				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The cigar has a start I stretch - checking if the database would have the same sequence: {}",&self.cigar);

				let mut cigar_vec = self.string_to_vec( &self.cigar );
				// this means likely that my read does reach into the < start area of the gene
				let count: usize = mat[1].parse().unwrap();

				for i in 0..count  { //  runs once for count==1
					if read.get_nucleotide_2bit( i ) == database.get_nucleotide_2bit( mapping_start - (count -i) ){
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} is the same as the sequence on read position {}", 
							mapping_start - (count -i), i );
						cigar_vec[i] = CigarEnum::Match;
					}else {
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} is NOT the same as the sequence on read position {}",
							 mapping_start - (count -i),
							 i
						);
						cigar_vec[i] = CigarEnum::Mismatch;
					}
				}
				ret = count;
				self.clear();
				self.reset_fom_path( &cigar_vec );
				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The updated Cigar looks like that: {self}");
			}
			// we now also need to check if the end would also contain I's
			let re_end = Regex::new(r"(\d+)I$").unwrap();
			if let Some(mat) =re_end.captures(&self.cigar) {

				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The cigar has a end I stretch - checking if the database would have the same sequence: {}",&self.cigar);

				let mut cigar_vec = self.string_to_vec( &self.cigar );
				// this means likely that my read does reach into the < start area of the gene
				let count: usize = mat[1].parse().unwrap();
				let (_mine, other) = self.calculate_covered_nucleotides(&self.cigar );
				let cigar_vec_len = cigar_vec.len();
				for i in 0..count  { //  runs once for count==1
					if read.get_nucleotide_2bit( read.len() - count + i ) == database.get_nucleotide_2bit( mapping_start + other +i  ) {
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} is the same as the sequence on read position {}",
							mapping_start + other +i  , read.len() - count + i +1
						);
						cigar_vec[ cigar_vec_len - i -1 ] = CigarEnum::Match;
					}else {
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} ({:?}) is NOT the same as the sequence on read position {} ({:?})", 
							mapping_start + other +i +1, database.get_nucleotide_2bit( mapping_start + other +i +1  ),
							read.len() - count + i , read.get_nucleotide_2bit( read.len() - count + i +1 )
						);
						cigar_vec[ cigar_vec_len - i -1 ] = CigarEnum::Mismatch;
					}
				}
				self.clear();
				self.reset_fom_path( &cigar_vec );
				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The updated Cigar looks like that: {self}");
			}
		}
		ret
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
				        let (mine, other) = self.calculate_covered_nucleotides(clipped_part);
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


	fn vec_to_cigar(&self, path: &[CigarEnum] ) -> String{
		let mut count = 0;
	    let mut last_direction = None;
	    let mut ret = "".to_string();
	    for &direction in path.iter() {
	        if Some(direction) == last_direction {
	            count += 1;
	            //println!("Count with state {last_direction:?} increased to {count}");
	        }
	        else {
	            if count > 0 {
	                ret += &format!("{}{}",count, last_direction.unwrap());
	            }
	            count = 1;
	            last_direction = Some(direction);
	        }
	    }
	    ret
	}

	/// converts a CigarEnum vector into a Cigar string and stores that internally.
	/// This function also updated the contains vector.
	pub fn reset_fom_path(&mut self, path: &[CigarEnum] ){
	    let mut count = 0;
	    let mut last_direction = None;

		self.cigar.clear();
		self.state_changes = 0;

		//println!("The cigar Path: {path:?}");

		//let mut loc_path = path.to_vec();

	    for &direction in path.iter() {
	        if Some(direction) == last_direction {
	            count += 1;
	            //println!("Count with state {last_direction:?} increased to {count}");
	        }
	        else {
	        	self.state_changes += 1;
	            if count > 0 {
	                self.cigar.push_str(&format!("{}{}",count, last_direction.unwrap()));
	                self.contains[last_direction.unwrap().to_id()] =  true ;
	                //println!("Added {count}{last_direction:?}");
	            }
	            count = 1;
	            last_direction = Some(direction);
	        }
	    }

	    if count > 0 {
	        self.cigar.push_str(&format!("{}{}",count, last_direction.unwrap()));
	        self.contains[last_direction.unwrap().to_id()] =  true ;
	    }
	}

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
	}



	pub fn calculate_cigar(&mut self, matrix : &Vec<Vec<Cell>>, last_match:bool )  {
			// Trace back the alignment path
		
		    let mut path = Vec::with_capacity( matrix.len().max( matrix[0].len() ));

		    let rows = matrix.len();
		    let cols =  matrix[0].len();

		    let mut i = rows - 1;
		    let mut j = cols - 1;

		    // fix the final entry
		    if  last_match {
				path.push(CigarEnum::Match);
		    }else {
		    	path.push(CigarEnum::Mismatch);
		    }

		    while i > 0 && j > 0 {
		        let current_score = matrix[i][j].score;
		        println!("Current score = {current_score}");
		        let diagonal_score = matrix[i - 1][j - 1].score;
		        let up_score = matrix[i - 1][j].score;
		        let left_score = matrix[i][j - 1].score;

		        let max = Self::max3(diagonal_score, up_score, left_score);
		        if max == diagonal_score {
		        	if diagonal_score > current_score {
		        		// mismatch!!!
		        		println!("adding Mismatch");
		        		path.push(CigarEnum::Mismatch);
		        	}else {
		        		// match
		        		println!("adding Match");
		        		path.push(CigarEnum::Match);
		        	}
		            i -= 1;
		            j -= 1;
		        } else if max == up_score {
		        	println!("adding Deletion");
		        	path.push(CigarEnum::Insertion);
		            i -= 1;
		        } else {
		        	println!("adding Insertion");
		            path.push(CigarEnum::Deletion);
		            j -= 1;
		        }
		    }

		    // If there are remaining gaps at the beginning of the sequences, fill them with the corresponding directions
		    while i > 0 {
		    	println!("End not reached by main - adding INSERTION");
		        path.push(CigarEnum::Insertion);
		        i -= 1;
		    }

		    while j > 0 {
		    	println!("End not reached by main - adding DELETION");
		        path.push(CigarEnum::Deletion);
		        j -= 1;
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