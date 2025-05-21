use std::fmt;
use regex::Regex;


#[derive(Clone, Copy, PartialEq, Debug)]
pub enum CigarEnum{
	Match,
	Mismatch,
	Insertion,
	Deletion,
	Empty,
	Softclip,
	Hardclip,
	Nothing,
}

impl CigarEnum{

	pub fn is_gap(&self) -> bool {
		self == &CigarEnum::Deletion || self == &CigarEnum::Insertion
	}

	pub fn is_clipp(&self) -> bool {
		self == &CigarEnum::Hardclip || self == &CigarEnum::Softclip
	}

	pub fn opposite(&self, other: &Self ) ->bool {
		other == &self.get_opposite()
	}

	/// Does this Cigar add to the databse length? 
	/// This can be used to calculate the length of the cigars in the end.
	pub fn adds_to_database(&self, report_n: bool ) -> bool {
		match &self{
			CigarEnum::Match => true,
			CigarEnum::Mismatch => true,
			CigarEnum::Insertion => false,
			CigarEnum::Deletion => true,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
			CigarEnum::Hardclip => false,
			CigarEnum::Softclip => false,
			CigarEnum::Nothing => report_n,
		}
	}

	/// Does this Cigar add to the read length? 
	/// This can be used to calculate the length of the cigars in the end
	pub fn adds_to_read(&self, report_n: bool ) -> bool {
		match &self{
			CigarEnum::Insertion => true,
			CigarEnum::Deletion => false,
			CigarEnum::Nothing => false,
			CigarEnum::Softclip => true,
			_ => self.adds_to_database( report_n )
		}
	}

	pub fn get_opposite(&self) -> Self {
		match &self{
			CigarEnum::Match => CigarEnum::Mismatch,
			CigarEnum::Mismatch => CigarEnum::Match,
			CigarEnum::Insertion => CigarEnum::Deletion,
			CigarEnum::Deletion => CigarEnum::Insertion,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
			CigarEnum::Hardclip => CigarEnum::Softclip,
			CigarEnum::Softclip => CigarEnum::Hardclip,
			CigarEnum::Nothing => panic!("You can not compare CigarEnum::Nothing to anything"),
		}
	}
	pub fn to_id(&self) -> usize{
		match &self{
			CigarEnum::Match => 0,
			CigarEnum::Mismatch => 1,
			CigarEnum::Insertion => 2,
			CigarEnum::Deletion => 3,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
			CigarEnum::Hardclip => 4,
			CigarEnum::Softclip => 5,
			CigarEnum::Nothing => 6,
		}
	}
	pub fn to_string(&self) -> String{
		match &self{
			CigarEnum::Insertion =>  "I".to_string(),
	        CigarEnum::Deletion => "D".to_string(),
	        CigarEnum::Match =>  "M".to_string(),
	        CigarEnum::Mismatch =>  "X".to_string(),
	        CigarEnum::Hardclip => "H".to_string(),
			CigarEnum::Softclip => "S".to_string(),
			CigarEnum::Nothing => "N".to_string(),
	        CigarEnum::Empty => panic!("That can not be convertet"),
	    }
	}
	pub fn from_str(c: &str) -> Option<CigarEnum> {
        match c {
            "I" => Some(CigarEnum::Insertion),
            "D" => Some(CigarEnum::Deletion),
            "M" => Some(CigarEnum::Match),
            "X" => Some(CigarEnum::Mismatch),
            "H" => Some(CigarEnum::Hardclip),
			"S" => Some(CigarEnum::Softclip),
			"N" => Some(CigarEnum::Nothing),
			// Add more cases as needed
            _ => None,
        }
    }
    pub fn get_regex() -> Regex {
    	Regex::new(r"(\d+)([MIDXHSN])").unwrap()
    }
}



// Implementing PartialEq for CigarEnum to allow comparison with a &str
impl PartialEq<&str> for CigarEnum {
    fn eq(&self, other: &&str) -> bool {
        match self {
            CigarEnum::Insertion => *other == "I",
            CigarEnum::Deletion => *other == "D",
            CigarEnum::Match => *other == "M",
            CigarEnum::Mismatch => *other == "X",
			CigarEnum::Hardclip => *other == "H",
			CigarEnum::Softclip => *other == "S",
			CigarEnum::Nothing => *other == "N",
            CigarEnum::Empty => false,  // We assume an empty variant should never compare to a string

        }
    }
}

impl fmt::Display for CigarEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string() )
    }
}