// /tests/analysis_gene_mapper_internals.rs

#[cfg(test)]
mod tests {

	use rustody::analysis_genomic_mapper::AnalysisGenomicMapper;
	use rustody::mapping_info::MappingInfo;
	use rustody::genes_mapper::SeqRec;
	use rustody::errors::MappingError;
	use rustody::singlecelldata::SingleCellData;

	use std::path::Path;
	use std::process::exit;


	fn test_this_seqence( seq: &[u8], database:String, sam_line: &str, err:Option<MappingError> ) -> (SingleCellData, Vec<String>){

		if !Path::new(&database).exists() {
	        eprintln!("The database file does not exist");
	        exit(1);
	    }

		let mut results = MappingInfo::new( None, 20.0, 10, None );
		let mut worker = AnalysisGenomicMapper::new( 32, "v1".to_string(), "mouse".to_string(),
			Some(database), 1, "bd", true);
		let pos = &[0,9, 21,30, 43,52, 52,60 ];

		worker.debug( Some(true) );
		// that contains a cell id for the version of the bd tool
		let r1 = SeqRec::new( b"SomeRead1", b"AGGAGATTAACTGGCCTGCGAGCCTGTTCAGGTAGCGGTGACGACTACATATGCTGCACATTTTTT", b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" );
	    // a read for the seq you wanted analyzed
	    let qual: Vec<u8>  = vec![b'F'; seq.len()];
	    let r2 = SeqRec::new( b"SomeRead2", seq, qual.as_slice() );
	    let data= vec![ (r1, r2 )];

	    let (single_cell_data, sam_strings) = worker.analyze_paralel( &data, &mut results, pos );
	    match err{
	    	Some(_e) => {
	    		assert_eq!( sam_strings.len(), 0, "I go no result for the search" );
	    		(single_cell_data, sam_strings)
	    	},
	    	None=> {
	    		if sam_strings.is_empty(){
	    			panic!("I go no result instead of a sam line!");
	    		}
	    		assert_eq!(sam_strings[0], sam_line, "We got the expected sam line?");
	    		(single_cell_data, sam_strings)
	    	}
	    }
	    
	}

    //A00681:881:H3MV7DSX7:1:1101:9679:1470 2:N:0:CACAATCCCA+TTGTGGATAT       0       Rps9    1       40      74M     *       0       0       TCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA     FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFF     NH:i:1  HI:i:1  AS:i:40 nM:i:0  RE:A:I  li:i:0  BC:Z:CGTACTTCGACA       QT:Z:FFFFFFFFFF:F       CR:Z:CCGATGGGTGAGAACC   CY:Z:FFFFFFFFFFFFF:FF   CB:Z:CCGATGGGTGAGAACC-1 UR:Z:CGTACTTCGACA       UZ:Z:FFFFFFFFFF:F       UB:Z:CGTACTTCGACA       RG:Z:Sample4:0:1:HN2CKBGX9:1

    //TCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA

    #[test]
    fn problematic_cigar() {
        let _name ="Rpl11_int";
        let database = "/home/med-sal/sens05_shared/common/genome/Rustody/version7/mouse/GRCm39_gm/".to_string();
        let _seq = b"---------------TCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA";
        let  _db = b"----ACAGGCAATGCTCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA";
		let  seq = b"GTTTGAAGGCAATGCTCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA";

        //let seq = b"TCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA";
        //let bam_line= "SomeRead2\t0\tRpl11_int\t320\t39\t1x84M\t*\t0\t0\tTCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
        let bam_line= "SomeRead2\t0\tRps9_int\t1\t39\t2X84M\t*\t0\t0\tTGAAGGCAATGCTCTCCTGCGGCGGCTTGTTCGCATTGGGGTGCTGGACGAGGGCAAGATGAAGCTGGATTACATCCTGGGCCTGAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.023255814\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
        let _ = test_this_seqence( seq, database, bam_line, None );
    }

}
