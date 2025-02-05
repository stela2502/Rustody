// test/creating_cigars.rs



#[cfg(test)]
mod tests {
    use rustody::analysis::AnalysisGeneMapper;
    use rustody::mapping_info::MappingInfo;
    use rustody::genes_mapper::SeqRec;
    use rustody::errors::MappingError;

    fn test_this_seqence( seq: &[u8], database:String, sam_line: Option<&str>, err:Option<MappingError> ){

        let mut results = MappingInfo::new( None, 20.0, 10, None );
        let mut worker = AnalysisGeneMapper::new( 32, "v1".to_string(), Some(database),
            None, "mouse".to_string(), None, 1, "bd", true);
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
            Some(e) => {
                assert_eq!( single_cell_data.is_empty(), true, "no results in the data object");
                assert_eq!( sam_strings.len(), 0, "I go no result for the search {:?}",e );
            },
            None=> {
                assert_eq!( single_cell_data.is_empty(), false, "there no result in the data object and I expected '{sam_line:?}'");
                match sam_line {
                    Some( sam ) => {
                        if sam_strings.is_empty(){
                            panic!("I go no result instead of a sam line!");
                        }
                        assert_eq!(sam_strings[0], sam, "We got the expected sam line?");
                    },
                    None => {
                        assert_eq!( sam_strings.len(), 0, "I go no result for the search" );
                    }
                }
                
            }
        }
        
    }

    #[test]
    fn test_deadlock1(){
        let seq = b"CCATTGCCCCCACGCTAGCTATATACTGAGGGAAGTGACCCTCCAGGGTTAGCTCAGATCTCTGATCGAACCCAC";
        let database = "testData/genes.fasta".to_string();

        let bam_line= "SomeRead2\t0\tLgals9\t0\t36\t12M1X31M1X5M1X17M2X3M1X1M\t*\t0\t0\tCCATTGCCCCCACGCTAGCTATATACTGAGGGAAGTGACCCTCCAGGGTTAGCTCAGATCTCTGATCGAACCCAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:36\tnM:i:0.08\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
        test_this_seqence( seq, database, Some(bam_line), None );
    }

    #[test]
    fn test_deadlock2(){
        let seq = b"ACAACATGGCTTCCAGAACAGTCGAGAGCAGAGTCTTGCCCCACCCACACCCATCCTGGAGGACAGTGGATAG";
        let database = "testData/genes.fasta".to_string();

        let bam_line= "SomeRead2\t0\tFam129c\t1\t37\t1S13M1D24M1X33M1X\t*\t0\t0\tACAACATGGCTTCCAGAACAGTCGAGAGCAGAGTCTTGCCCCACCCACACCCATCCTGGAGGACAGTGGATAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:37\tnM:i:0.054054055\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
        test_this_seqence( seq, database, Some(bam_line), None );
    }

}