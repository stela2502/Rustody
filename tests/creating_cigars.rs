// test/creating_cigars.rs



#[cfg(test)]
mod tests {
	//use rustody::genes_mapper::Cigar;
	use rustody::genes_mapper::gene_data::GeneData;
	//use rustody::traits::BinaryMatcher;
	use rustody::genes_mapper::NeedlemanWunschAffine;

	#[test]
	fn test_full_match(){
		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let seq2 = &seq1.clone();
		let gd1 = GeneData::new( seq1, "read", "read", "chrM", 0, false );
		let gd2 = GeneData::new( seq2, "database", "database", "chrM", 0,true);

		let mut nwa = NeedlemanWunschAffine::new();
		nwa.set_debug(true);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4);

		assert_eq!( &nwa.cigar.cigar, "90M", "I expected 90M as the sequences are the same and they are 90 bp long" );
	}

	#[test]
	fn test_deletion_match(){

		//
		//
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATA----------------------------------TTATTTTACTTCTACTCACAAT
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT

		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATATTATTTTACTTCTACTCACAAT";
		let seq2 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let gd1 = GeneData::new( seq1, "read", "read", "chrM", 0 , false );
		let gd2 = GeneData::new( seq2, "database", "database", "chrM", 0, true );

		let mut nwa = NeedlemanWunschAffine::new();
		nwa.set_debug(true);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4);

		//assert_eq!( &cigar.to_string(), "67D23M", "I expected 34M34D22M as I manually deleted 34 bp from the read" );
		assert_eq!( &nwa.cigar.cigar, "32M34D24M", "I expected 33M37D20M as I manually deleted 34 bp from the read\n{}", nwa.cigar.as_alignement( &gd1, &gd2) );
	}

	#[test]
	fn test_insertion_match(){
		//
		//
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATA--------------------------------- TTATTTTACTTCTACTCACAAT

		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let seq2 =                                   b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATATTATTTTACTTCTACTCACAAT";

		let gd1 = GeneData::new( seq1, "read", "read", "chrM", 0, false );
		let gd2 = GeneData::new( seq2, "database", "database", "chrM", 0, true );

		let mut nwa = NeedlemanWunschAffine::new();
		nwa.set_debug(true);
		
		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4 );
	
		//assert_eq!( &cigar.to_string(), "67I23M", "I expected 34M34I22M as I manually deleted 34 bp from the database" );
		assert_eq!( &nwa.cigar.cigar, "32M34I24M", "I expected 32M34I24M as I manually deleted 34 bp from the database" );
	}
}