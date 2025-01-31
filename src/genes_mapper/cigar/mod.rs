// genes_mapper/cigar/mod.rs

pub mod cigar;
pub mod cigar_enum;
pub mod cigar_tuple;

pub use cigar::Cigar as Cigar;
pub use cigar::CigarEndFix as CigarEndFix;
pub use cigar_enum::CigarEnum as CigarEnum;

pub use cigar_tuple::CigarTuple as CigarTuple;