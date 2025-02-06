// /tests/cigar.rs


#[cfg(test)]
mod tests {
	use rustody::genes_mapper::Cigar;
	use rustody::genes_mapper::cigar::{CigarEnum, CigarEndFix, CigarTuple};
	//easiest class with BinaryMatcher
	use rustody::genes_mapper::GeneData;
    use rustody::traits::BinaryMatcher;

    #[test]
    fn simple_stupid_clone_test(){
    	let mut original = Cigar::new( "30M" );
    	original.dropped_start = 5;
    	original.dropped_end = 10;
    	let copy = original.clone();
    	assert_eq!( original, copy, "Clone copied also the dropped values");
    }

	#[test]
	fn test_is_good_gap() {
		let a = GeneData::new( b"AAACTGTTT", "a", "a_", "chr1", 1, false);
		let b = GeneData::new( b"AAACGTTT", "b", "b_", "chr1", 1, false);
		let cigar = Cigar::new( "4M1I4M" );
		let vec = cigar.to_cigar_tupel_vec(false);
		assert_eq!( cigar.is_good_gap( &vec[1], &a, &b ), true, "Insert gap is good? {}", cigar.as_alignement(&a, &b) );
		let c = GeneData::new( b"TTTCGAAA", "c", "c_", "chr1", 1, false);
		assert_eq!( cigar.is_good_gap( &vec[1], &a, &c ), false, "Insert gap is bad? {}", cigar.as_alignement(&a, &c) );

		let cigar2 = Cigar::new( "4M1D4M" );
		let vec2 = cigar2.to_cigar_tupel_vec(false);
		assert_eq!( cigar2.is_good_gap( &vec2[1], &b, &a ), true, "Deletion gap is good? {}", cigar2.as_alignement(&b, &a) );
		assert_eq!( cigar2.is_good_gap( &vec2[1], &c, &a ), false, "Deletion gap is bad? {}", cigar2.as_alignement(&c, &a) );
	}

	#[test]
	fn test_is_bad_long_gap() {
		let a = GeneData::new( b"TGGGGCCCCTGTCACTC", "a", "a_", "chr1", 1, false);
	    let b = GeneData::new( b"TGGGGTCCCTGTCACTC", "b", "b_", "chr1", 1, false);
	    let cigar = Cigar::new( "5M4D3I1M1I7M" );
	    let vec = cigar.to_cigar_tupel_vec(false);

	    assert_eq!( cigar.is_good_gap( &vec[1], &a, &b ), false, "Insert gap is bad? {}", cigar.as_alignement(&a, &b) );
	}


	#[test]
	fn test_change_fixed_manually(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("5D30M");
		obj.fixed = Some( CigarEndFix::Both );
		assert_eq!( obj.fixed, Some( CigarEndFix::Both ), "Manual chjange works");
	}

	#[test]
	fn test_cvalculate_covered_1(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("5D30M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (30, 35), "corect sizes" );

		obj.restart_from_cigar("5I30M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (35, 30), "corect sizes" );

		obj.restart_from_cigar("5I15M3D15M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (35, 33), "corect sizes" );

		obj.restart_from_cigar("5D15M3I15M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (33, 35), "corect sizes" );
	}

	#[test]
	fn test_to_sam_string(){
		let mut obj = Cigar::default();

		obj.restart_from_cigar("1D15M1I30M");
		let fixed = obj.to_sam_string(46);

		assert_eq!( fixed, ("1D15M1I30M".to_string(),0 ), "internal insert overhanging D fixed at start");
		
		obj.restart_from_cigar("15M1I30M1D");
		let fixed = obj.to_sam_string(46);
		assert_eq!( fixed, ("15M1I30M1D".to_string(),0 ), "internal insert overhanging D fixed at end");
	}

	#[test]
	fn test_soft_start(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M");
		obj.soft_clip_start_end( );
		println!("This is the obtained cigar: {obj}");
		assert_eq!( obj.cigar, "24X65M");
		assert_eq!( obj.fixed, Some(CigarEndFix::Start), "start fixed");
	}

	#[test]
	fn test_soft_end(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("59M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		obj.soft_clip_start_end( );
		println!("This is the obtained cigar: {obj}");
		assert_eq!( obj.cigar, "59M30X");
		assert_eq!( obj.fixed, Some(CigarEndFix::End), "end fixed");
	}

	#[test]
	fn test_soft_both(){
		let mut obj = Cigar::default();
		println!("default cigar: {:?}", obj);
		obj.restart_from_cigar("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		println!("restart_from_cigar: {:?}", obj);
		println!("Memory address of obj before function: {:p}", &obj);
		obj.soft_clip_start_end();
		println!("Immediate post-function fixed: {:?}", obj.fixed);
		println!("Memory address of obj after function: {:p}", &obj);
		println!("This is the obtained cigar: {obj}");
		assert_eq!( obj.cigar, "24X65M30X");
		assert_eq!( obj.fixed, Some(CigarEndFix::Both), "both fixed");
	}

	#[test]
	fn test_quality(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		obj.soft_clip_start_end();
		assert_eq!( obj.mapping_quality(), 21 );
	}
	#[test]
	fn test_cigar_fix(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1M1X1M5X1M1X2M2X1M1X1M1X2M8X2I61M");
		obj.soft_clip_start_end();
		println!("This is the obtained cigar: {obj}");
		assert_eq!( obj.cigar, "30X61M");
		assert_eq!( obj.fixed, Some(CigarEndFix::Start), "both fixed");
	}


	#[test]
	fn test_calculate_covered_nucleotides() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("14M1X39M1X7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (75, 75), "corect sizes" )
	}

	#[test]
	fn test_calculate_covered_nucleotides_deletions() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("14M1X39M1D7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (74, 75), "corect sizes" )
	}
	#[test]
	fn test_calculate_covered_nucleotides_insertions() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("14M1X39M1I7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (75, 74), "corect sizes" )
	}
	#[test]
	fn test_calculate_covered_nucleotides_real() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X8M1I39M2X16M7X");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (74, 73), "corect sizes" )
	}

	#[test]
	fn test_calculate_len_real() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X8M1I39M2X16M7X");
		assert_eq!( obj.len(), 74, "1X8M1I39M2X16M7X - corect sizes" )
	}


	#[test]
	fn test_default_is_worst() {
		let mut obj1 = Cigar::default();
		let obj2 = Cigar::default();

		obj1.restart_from_cigar( "32M" );

		assert!( obj1.better_as(&obj2), "{obj1:?}\nis better than \n{obj2:?} ({})", obj1.better_as(&obj2) );
	}

	#[test]
	fn test_compare() {
		let mut obj1 = Cigar::new("32M");
		let mut obj2 = Cigar::new("36M");
		//obj1.reset_fom_path(&vec![CigarTuple::from_scratch( CigarEnum::Match, 32)]);
		//obj2.reset_fom_path(&vec![CigarTuple::from_scratch( CigarEnum::Match, 32)]);

		assert!( obj2.better_as(&obj1), "#1 \n{obj2} is better than \n{obj1}? {}",  obj2.better_as(&obj1) );

		obj1.reset_fom_path( &[
			CigarTuple::from_scratch( CigarEnum::Insertion, 1), 
			CigarTuple::from_scratch( CigarEnum::Match,6), 
			CigarTuple::from_scratch( CigarEnum::Deletion, 1), 
			CigarTuple::from_scratch( CigarEnum::Match, 4) 
			] 
		);
		obj2.reset_fom_path( &[
			CigarTuple::from_scratch( CigarEnum::Match,2),
			CigarTuple::from_scratch( CigarEnum::Deletion, 1),
			CigarTuple::from_scratch( CigarEnum::Match,2),
			CigarTuple::from_scratch( CigarEnum::Deletion, 1),
			CigarTuple::from_scratch( CigarEnum::Match,1),
			CigarTuple::from_scratch( CigarEnum::Insertion,1),
			CigarTuple::from_scratch( CigarEnum::Match,4),
		    ]
		);
		assert!( obj1.better_as(&obj2), "#2 \n{obj1:?} is better than \n{obj2:?} ({})", obj2.better_as(&obj1) );

		obj1.restart_from_cigar( "21M1D49M1I3M" ); // len 74
		obj2.restart_from_cigar( "1I20M1D52M" ); // len 73

		assert!( obj1.better_as(&obj2), "#3 \n{obj1:?}\nis better than \n{obj2:?} ({})", obj1.better_as(&obj2) );
	}

	



}