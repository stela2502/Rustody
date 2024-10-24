// genes_mapper/mod.rs

pub mod genes_mapper;
pub mod gene_data;
pub mod gene_link;
pub mod mapper_result;
pub mod sequence_record;
pub mod multi_match;
pub mod needleman_wunsch_affine;

pub mod cigar;

pub use genes_mapper::GenesMapper as GenesMapper;
pub use cigar::Cigar as Cigar;
pub use crate::genes_mapper::mapper_result::MapperResult as MapperResult;
pub use crate::genes_mapper::sequence_record::SeqRec as SeqRec;
pub use needleman_wunsch_affine::NeedlemanWunschAffine as NeedlemanWunschAffine;

pub use crate::genes_mapper::multi_match::MultiMatch as MultiMatch;
pub use cigar::CigarEndFix as CigarEndFix;
