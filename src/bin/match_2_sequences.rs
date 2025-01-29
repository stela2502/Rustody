use clap::Parser;

use rustody::genes_mapper::gene_data::GeneData;
use rustody::genes_mapper::needleman_wunsch_affine::NeedlemanWunschAffine;
use rustody::traits::BinaryMatcher;


#[derive(Parser)]
#[clap(version = "0.1.1", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    read: String,
    /// the input R2 samples file
    #[clap(short, long)]
    database: String,
}


fn main() {
    
    let opts: Opts = Opts::parse();
    // unique_name:&str, name:&str, chr:&str, start:usize, index_type: bool
    let read = GeneData::new( opts.read.as_bytes(), "read", "read", "chr1", 1, false );
    let database =  GeneData::new( opts.database.as_bytes(), "database", "database", "chr21", 1, false );

    let mut matcher = NeedlemanWunschAffine::new();
    matcher.initialize( database.len(), read.len());

    matcher.needleman_wunsch_affine( &read, &database, 10.0 );

    println!( "{}",matcher.to_string(&read, &database, 10.0 ));

}
