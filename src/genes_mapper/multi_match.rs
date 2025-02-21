use crate::genes_mapper::MapperResult;
//use crate::genes_mapper::{ CigarEndFix};

//use rand::Rng;

use core::fmt;

pub struct MultiMatch{
	data:Vec<MapperResult>
}


// Implementing Display trait for MultiMatch
impl fmt::Display for MultiMatch {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		let formatted_data: String = self.data.iter().map(|x| format!("{:?}, len: {}", x, x.get_cigar().len())).collect::<Vec<_>>().join("\n");
		write!(f, "MultiMatch has {} entries: ( \n{} )", self.data.len(), formatted_data )
	}
}



impl MultiMatch{
	pub fn new() -> Self{
		Self{
			data: Vec::<MapperResult>::with_capacity(4),
		}
	}
	pub fn push( &mut self, data:MapperResult) {
		self.data.push( data );
	}
	pub fn clear( &mut self ) {
		self.data.clear();
	}
	pub fn len( &self ) -> usize{
		self.data.len()
	}

	pub fn get_best( &mut self, _length:usize ) -> Result<MapperResult, MapperResult > {
		//println!( "MultiMatch::get_best has been called! {self}");
		self.data.sort_unstable_by(|a, b| b.cigar().cmp(&a.cigar())); // Sort in place

		// println!("My sorted results - best first?:\n{:?}",self.data);
	    match self.data.as_slice() {
	        [best, second, ..] if best == second => Err(best.clone()), // Tie case now also incluging the gene name differences (sorting _int after!)
	        [best, ..] => Ok(best.clone()), // Unique best
	        _ => Err(MapperResult::default()), // No results case (handle appropriately)
	    }
	}

}