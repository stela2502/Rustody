// mapper_result.rs

use crate::genes_mapper::Cigar;
//use crate::genes_mapper::CigarEndFix;
//use crate::genes_mapper::gene_data::GeneData;

use std::cmp::Ordering;

use core::fmt;

#[derive(Debug, Clone )]
pub struct MapperResult{
	gene_id: usize,
	start: usize,
	multimapper: bool,
	cigar: Option<Cigar>,
	mapq: u8,
	score: usize,
    nw: f32,
	edit_dist: f32,
	gene_name:String,
    /// the length of the database entry - is necessary to estimate if a short match is good enough (e.g. sampleid or AB tag)
    db_length: usize,
    // has the read been clipped to obtain this match?!
    //clipped: (usize, usize),

}

// Implementing Display trait for MapperResult
impl fmt::Display for MapperResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MapperResult ( gene_id {}, gene_name {}, start {}, multimapper {}, cigar {:?}, nw {}, mapq {}, edit_dist {}, position_from_end {})\n", 
            self.gene_id, self.gene_name, self.start, self.multimapper, self.cigar, self.nw, self.mapq, self.edit_dist, self.position_from_end() )
    }
}

impl Default for MapperResult {
    fn default() -> Self {
        MapperResult {
            gene_id: 0,
            start: 0,
            multimapper: false,
            cigar: Some(Cigar::default()),
            mapq: 200,
            score: 0,
            nw:100.0,
            edit_dist: 100.0,
            gene_name: String::default(),
            db_length: 0,
        }
    }
}


impl PartialEq for MapperResult {
    fn eq(&self, other: &Self) -> bool {
        self.cigar == other.cigar && self.gene_name == other.gene_name
    }
}

impl Eq for MapperResult {}

impl PartialOrd for MapperResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other)) // Delegate to `Ord`
    }
}

/// this will automaticall take care of the same cigar, but attached _int in the gene names.
impl Ord for MapperResult {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.cigar.cmp(&other.cigar) {
            Ordering::Equal => {        
                // Otherwise, use cigar comparison 
                self.gene_name.cmp(&other.gene_name)// Compare gene_name if cigars are equal
            },
            other => other,
        }
    }
}


impl MapperResult{
	pub fn new( gene_id:usize, start: usize, multimapper:bool, cigar:Option<Cigar>,
        mapq:u8, nw:f32,score:usize, edit_dist:f32, name:&str, db_length:usize, 
     ) -> Self{

        Self{
			gene_id,
			start,
			multimapper,
			cigar,
			mapq,
			score,
            nw,
			edit_dist,
			gene_name: name.to_string(),
            db_length,
		}
	}

    pub fn db_length(&self) -> usize{
        self.db_length
    }

    pub fn get_nw(&self) -> f32{
        self.nw
    }
	pub fn get_name(&self) -> &str{
		&self.gene_name
	}

	pub fn gene_id(&self) -> usize {
        self.gene_id
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn multimapper(&self) -> bool {
        self.multimapper
    }

    pub fn cigar(&self) -> &Option<Cigar> {
        &self.cigar
    }

    pub fn get_cigar(&self) -> Cigar{
        match &self.cigar{
            Some(cig) => cig.clone(),
            None => Cigar::default(),
        }
    }

    pub fn position_from_end (&self) -> usize{
        self.db_length.saturating_sub( self.start )
    }

    pub fn mapq(&self) -> u8 {
        self.mapq
    }
    pub fn score(&self) -> usize {
        self.score
    }
    pub fn edit_dist(&self) -> f32 {
        self.edit_dist
    }

    pub fn save(&self) -> bool {
        true
    }
}