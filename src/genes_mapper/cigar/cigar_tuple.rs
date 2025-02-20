
use crate::traits::BinaryMatcher;
use crate::genes_mapper::cigar::CigarEnum;


use std::fmt;


#[derive(Debug, Clone)]
pub struct CigarTuple {
	/// the start position in the CigarEnum Vector
	pub vec_pos: usize,
	/// the start position in the Cigar string
	pub str_pos:usize,
	/// the length in CigarEnum Vector dimension
	pub vec_len : usize,
	/// the length in CigarString dimension
	pub str_len:usize,
	/// the position on the read at this position
	pub read_position:usize,
	/// the position on the database at this position
	pub database_position:usize,
	/// The Cigar type of this operation
	pub option:CigarEnum,
}


impl CigarTuple {
    /// Creates a CigarTuple from a regex capture group.
    pub fn from_match(operation: &str, num_str: &str, start_on_vec: usize, start_in_str: usize, read_position: usize, database_position:usize ) -> Option<Self> {
        // Capture the length and operation from the capture groups

        // Parse the length from the captured string
        let length: usize = num_str.parse::<usize>().ok()?;

        // Convert the operation string to a CigarEnum
        let cigar_type = CigarEnum::from_str(operation).unwrap_or_else(|| {
            panic!(
                "Invalid CIGAR operation '{}' in the captured string.",
                operation
            )
        });

        // Return a new CigarTuple with the parsed data
        Some(Self {
            vec_pos: start_on_vec,
            str_pos: start_in_str,
            vec_len: length,
            str_len: num_str.len(),
            read_position,
            database_position,
            option: cigar_type,
        })
    }

    pub fn from_scratch( cigar:CigarEnum, length:usize ) ->Self {
    	Self {
            vec_pos: 0,
            str_pos: 0,
            vec_len: length,
            str_len:  ((length as f64).log10().floor() as usize) + 2,
            read_position:0,
            database_position:0,
            option: cigar,
        }
    }
 	pub fn print_debug(&self) {
    	println!("{} vec {}..{}; str {}..{} on_read: {} on_database:{}",
    		self,
    		self.vec_pos, self.vec_pos+ self.vec_len,
    		self.str_pos, self.str_pos + self.str_len,
    		self.read_position,
    		self.database_position
    		)
    }

    pub fn slice_from_database<T>(&self, database:&T ) -> Option<String>where
    T: BinaryMatcher{
    	if self.option == CigarEnum::Insertion{
    		database.as_str( self.read_position, self.vec_len )
    	}else {
    		database.as_str( self.database_position, self.vec_len )
    	}
    	
    }
    pub fn slice_from_read<T> (&self, read:&T) -> Option<String> where
    T: BinaryMatcher{
    	if self.option == CigarEnum::Deletion{
    		read.as_str( self.database_position, self.vec_len )
    	}else {
    		read.as_str( self.read_position, self.vec_len )
    	}
    }
    pub fn str_start(&self) -> usize{
    	self.str_pos
    }
    pub fn str_end( &self ) -> usize{
    	self.str_pos + self.str_len
    }
    pub fn vec_start(&self) -> usize{
    	self.vec_pos
    }
    pub fn vec_end(&self) -> usize{
    	self.vec_pos + self.vec_len
    }
    pub fn is_gap(&self) ->bool {
    	self.option.is_gap()
    }
    // this will return the vector lengh as this is what is intereting most of the times anyhow.
    pub fn len(&self) -> usize {
    	self.vec_len
    }

    pub fn extend_vec(&self, vec: &mut Vec<CigarEnum>) {
        vec.extend(std::iter::repeat(self.option).take(self.len()));
    }

    /// This will 'shorten' the own and return the one that should be inserted after this one
    /// as a new object.
    pub fn split_at (&mut self, rel_pos:usize, new_option: CigarEnum ) -> Option<Self>{
        if self.len() < rel_pos {
            self.option = new_option;
            None
        }else {
            self.vec_len -= rel_pos;
            Some(Self::from_scratch(new_option, rel_pos ))
        }
    }

    pub fn to_string(&self) -> String{
        format!("{}{}", self.vec_len, self.option )
    }

    pub fn is_a(&self, other:&CigarEnum ) -> bool {
        &self.option == other
    }


    // You can add other methods to operate on the `CigarTuple`

}


impl fmt::Display for CigarTuple {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    	if self.len() == 0{
    		 write!(f, "")
    	}else {
    		write!(f, "{}{}", self.vec_len, self.option )
    	}
    }
}
