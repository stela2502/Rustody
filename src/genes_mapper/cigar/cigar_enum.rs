use std::fmt;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum CigarEnum{
	Match,
	Mismatch,
	Insertion,
	Deletion,
	Empty,
}

impl CigarEnum{

	pub fn is_gap(&self) -> bool {
		self == &CigarEnum::Deletion || self == &CigarEnum::Insertion
	}
	pub fn opposite(&self, other: &Self ) ->bool {
		match &self{
			CigarEnum::Match => other == &CigarEnum::Mismatch,
			CigarEnum::Mismatch => other == &CigarEnum::Match,
			CigarEnum::Insertion => other == &CigarEnum::Deletion,
			CigarEnum::Deletion => other == &CigarEnum::Insertion,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
		}
	}
	pub fn get_opposite(&self) -> Self {
		match &self{
			CigarEnum::Match => CigarEnum::Mismatch,
			CigarEnum::Mismatch => CigarEnum::Match,
			CigarEnum::Insertion => CigarEnum::Deletion,
			CigarEnum::Deletion => CigarEnum::Insertion,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
		}
	}
	pub fn to_id(&self) -> usize{
		match &self{
			CigarEnum::Match => 0,
			CigarEnum::Mismatch => 1,
			CigarEnum::Insertion => 2,
			CigarEnum::Deletion => 3,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
		}
	}
	pub fn to_string(&self) -> String{
		match &self{
			CigarEnum::Insertion =>  "I".to_string(),
	        CigarEnum::Deletion => "D".to_string(),
	        CigarEnum::Match =>  "M".to_string(),
	        CigarEnum::Mismatch =>  "X".to_string(),
	        CigarEnum::Empty => panic!("That can not be convertet"),
	    }
	}
	pub fn from_str(c: &str) -> Option<CigarEnum> {
        match c {
            "I" => Some(CigarEnum::Insertion),
            "D" => Some(CigarEnum::Deletion),
            "M" => Some(CigarEnum::Match),
            "X" => Some(CigarEnum::Mismatch),
            // Add more cases as needed
            _ => None,
        }
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
            CigarEnum::Empty => false,  // We assume an empty variant should never compare to a string
        }
    }
}

impl fmt::Display for CigarEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let direction_str = match self {
            CigarEnum::Match => "M",
			CigarEnum::Mismatch => "X",
			CigarEnum::Insertion => "I",
			CigarEnum::Deletion => "D",
			CigarEnum::Empty => panic!("There is an empty cigar entry in your vector!"),
        };
        write!(f, "{}", direction_str)
    }
}