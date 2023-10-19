use clap::Parser;

use regex::Regex;

use this::fast_mapper::FastMapper;
use this::gene_family::GeneFamily;
use this::gene::Gene;

use needletail::parse_fastx_file;

//use this::sampleids::SampleIds;
//use this::analysis::

//use std::path::PathBuf;
use std::fs;
//use std::path::Path;

use std::time::SystemTime;

use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::collections::HashSet;

// use std::convert::TryInto;

/// Just run a test case for speed reproducability and simplicity

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the gff gene information table
    #[clap(default_value= "testData/TE_quant/chr14_KI270722v1_random_TEtranscripts.gtf.gz",short, long)]
    gtf: String,
    /// the fasta genome data
    #[clap( default_value=  "testData/TE_quant/chr14_KI270722v1_random.fa.gz",short, long)]
    file: String,
    /// the outpath
    #[clap(default_value=  "testData/TEmapperTest",short, long)]
    outpath: String,
    #[clap(default_value= "testData/MyAbSeqPanel.fasta", short, long)]
    antibody: String,
    /// the mapping kmer length
    #[clap(default_value_t=32, long)]
    gene_kmers: usize,
    /// create text outfile instead of binary
    #[clap(default_value_t=false, long)]
    text: bool,
}

/*
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/

// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    let mut kmer_size = opts.gene_kmers;
    if opts.gene_kmers > 32{
        eprintln!("Sorry the max size of the kmers is 32 bp");
        kmer_size = 32;
    }

    const COVERED_AREA:usize = 400; // cover 600 bp of the transcript
    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

    // let log_file_str = PathBuf::from(&opts.outpath).join(
    //     "index_log.txt"
    // );

    // println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );
    
    // let log_file = match File::create( log_file_str ){
    //     Ok(file) => file,
    //     Err(err) => {
    //         panic!("Error: {err:#?}" );
    //     }
    // };


    // read the fasta data in!
    let mut seq_records = HashMap::< String, Vec<u8>>::new();

    let mut expr_file = parse_fastx_file(&opts.file).expect("valid path/file");
    let mut id: &[u8];
    let delimiter: &[u8] = b"  ";  // The byte sequence you want to split by

    while let Some(e_record) = expr_file.next() {
        let seqrec = e_record.expect("invalid record");
        id = seqrec.id().split(|&x| x == delimiter[0]).collect::<Vec<&[u8]>>()[0];
        //id = seqrec.id().split(|&x| x == delimiter[0]).collect()[0];
        eprintln!("'{}'", std::str::from_utf8(id).unwrap() );
        seq_records.insert( std::str::from_utf8( id ).unwrap().to_string() , seqrec.seq().to_vec());
    }




    // and init the Index:

    let mut index = FastMapper::new( kmer_size );

    // in short I need to get an internal model of a gene to work.
    // I want to know where the gene starts ans ends (likely transcripts)
    // I want to get intron/exon bounds and then add both a spliced version of the transcript 
    // as well as an unspliced version

    /*
    seqname   0
    source    1    
    Feature   2
    start     3
    end       4
    score     5
    strand    6
    frame     7
    attribute 8 
    */
    
    let mut genes = HashMap::<String, GeneFamily>::new();

    let gtf = Regex::new(r".*gtf.?g?z?$").unwrap();
    let chr = Regex::new(r"^chr").unwrap();

    let re_gene_name: Regex;
    let re_gene_id: Regex;
    let re_transcript_id : Regex;


    match gtf.is_match( &opts.gtf ){
        true => {
            eprintln!("gtf mode");
            re_gene_name =  Regex::new(r#".* gene_name "([\(\)\w\d\-\._]*)""#).unwrap();
            re_gene_id = Regex::new(r#"gene_id "([\(\)\w\d\-\._]*)";"#).unwrap();
            re_transcript_id = Regex::new(r#"transcript_id "([\(\)\w\d\-\._]*)";"#).unwrap();
        },
        false => {
            eprintln!("gff mode.");
            re_gene_name =  Regex::new(r#".* ?gene_name=([\(\)\w\d\-\._]*); ?"#).unwrap();
            re_gene_id = Regex::new(r#"gene_id=([\(\)\w\d\-\._]*);"#).unwrap();
            re_transcript_id = Regex::new(r#"transcript_id=([\(\)\w\d\-\._]*);"#).unwrap();
        },
    }
    //let mut gene_id:String;
    let mut gene_name:String;
    let mut transcript_id:String;

    // now we need to read the gtf info and get all genes out of that
    // we specificly need the last 64bp of the end of the gene as that is what we are going to map
    // if there is a splice entry in that area we need to create a spliced and an unspliced entry

    if ! opts.gtf.ends_with(".gz") {
        panic!("Please gzip your gtf file - thank you! {}", opts.gtf.to_string());
    }

    let f1 = match File::open( &opts.gtf ){
        Ok(file) => file,
        Err(err) => panic!("The file {} does not exists: {err}", &opts.gtf ),
    };
    let file1 = GzDecoder::new(f1);
    let reader = BufReader::new( file1 );
    let mut missing_chr:HashSet<String>  = HashSet::new();

    for line in reader.lines() {
        let rec = line.ok().expect("Error reading record.");
        let mut parts: Vec<String> = rec.split('\t').map(|s| s.to_string()).collect();
        if parts.len() < 8{
            continue;
        }
        // if parts[2] == "transcript"{
        //     // capture the parts I need using my regexp modules
        //     if let Some(captures) = re_gene_name.captures( &parts[8].to_string() ){
        //         gene_name = captures.get(1).unwrap().as_str().to_string();
        //     }else {
        //         if let Some(_captures) = re_gene_id.captures( &parts[8].to_string() ){
        //             continue; // this likely clutters up the data anyhow.
        //         }
        //         else {
        //             panic!("I could not identify a gene_name in the attributes {:?}", &parts[8].to_string() );
        //         }
        //     }
        //     // if let Some(captures) = re_gene_id.captures( &parts[8].to_string() ){
        //     //     gene_id = captures.get(1).unwrap().as_str().to_string();
        //     // }else {
        //     //     panic!("I could not identify a gene_id in the attributes {:?}", &parts[8].to_string() );
        //     // }
        //     if let Some(captures) = re_transcript_id.captures( &parts[8].to_string() ){
        //         transcript_id = captures.get(1).unwrap().as_str().to_string();
        //     }else {
        //         panic!("I could not identify a transcript_id in the attributes {:?}", &parts[8].to_string() );
        //     }
            
        //     // and add a gene
        //     // pub fn new(chrom:String, start_s:String, end_s:String, sense_strand_s:String, name:String, id:String )
        //     let gene = Gene::new( parts[0].to_string(),  parts[3].to_string(), parts[4].to_string(), parts[6].to_string(), gene_name.to_string(), transcript_id.to_string() );
        //     genes.insert( transcript_id.clone(), gene );
        // }

        if parts[2] == "exon"{
            // capture the parts I need
            if let Some(captures) = re_transcript_id.captures( &parts[8].to_string() ){
                transcript_id = captures.get(1).unwrap().as_str().to_string();
                // and add an exon
                match genes.get_mut( &transcript_id.clone().to_string() ){
                    Some(gene) => gene.add_exon( parts[3].to_string(), parts[4].to_string() ),
                    None => eprintln!( "ignoring transcript! ({})", transcript_id.to_string())
                }
                // Ha this is a brutal hack for the TE gtfs only - just treat this entry as a gene, too.
                if let Some(captures) = re_gene_name.captures( &parts[8].to_string() ){
                    gene_name = captures.get(1).unwrap().as_str().to_string();
                    let gene = Gene::new( parts[0].to_string(),  parts[3].to_string(), parts[4].to_string(), parts[6].to_string(), gene_name.to_string(), transcript_id.to_string() );
                    genes.insert( transcript_id, gene );
                }else {
                    panic!("I could not identify a gene_name in the attributes {:?}", &parts[8].to_string() );
                }
                // and issue the clean up
                parts[2] = "gene".to_string();
            }
            
        }

        if parts[2] == "gene"{
            // here we should start to 'clean out the old ones'!
            let start = match parts[3].parse::<usize>(){
                Ok(v) => v,
                Err(e) => panic!("I could not parse the start of the transcript as usize: {e:?}"),
            };
            let mut to_remove = Vec::new();

            // clean up
            'genes: for (k, gene) in &genes {
                if let Some(_entry) = missing_chr.get( &gene.chrom.to_string() ){
                    eprintln!("gene '{}' - no sequence for chr '{}'", &gene.name.to_string(), &gene.chrom.to_string() );
                    continue 'genes;
                }
                if gene.passed(start) {
                    // Do something with the gene, e.g. remove it
                    match seq_records.get( &gene.chrom.to_string() ){
                        Some(seq) => {
                            gene.add_to_index( seq, &mut index, COVERED_AREA );
                            //println!("The genes detected: {:?}", index.names_store );
                        },
                        None => {
                            if chr.is_match ( &gene.chrom.to_string() ){
                                match seq_records.get( &gene.chrom.to_string()[3..] ){
                                    Some(seq) => {
                                        gene.add_to_index( seq, &mut index, COVERED_AREA );
                                        //println!("The genes detected: {:?}", index.names_store );
                                    },
                                    None => {
                                        missing_chr.insert( gene.chrom.to_string() );
                                        eprintln!("I do not have the sequence for the chromosome {}", gene.chrom.to_string() );
                                    }
                                }
                            }else {
                                match seq_records.get( &format!("chr{}", &gene.chrom.to_string()) ){
                                    Some(seq) => {
                                        gene.add_to_index( seq, &mut index, COVERED_AREA );
                                        //println!("The genes detected: {:?}", index.names_store );
                                    },
                                    None => {
                                        missing_chr.insert( gene.chrom.to_string() );
                                        eprintln!("I do not have the sequence for the chromosome {}", gene.chrom.to_string() );
                                    }
                                }

                            }
                        }
                            
                    }
                    
                    to_remove.push(k.clone());
                }
            }

            for k in to_remove {
                genes.remove(&k);
            }
        }
    }

    for (_, gene) in &genes {
        // Do something with the gene, e.g. remove it
        match seq_records.get( &gene.chrom.to_string() ){
            Some(seq) => {
                gene.add_to_index( seq, &mut index, COVERED_AREA );
                //println!("The genes detected: {:?}", index.names_store );
            },
            None => {
                let keys_vec: Vec<_> = seq_records.keys().collect();
                panic!("I do not have the sequence for the chromosome {} {:?}", gene.chrom.to_string(), keys_vec )
            },
        }

    }
    eprintln!(" total first keys {}\n total second keys {}\n total single gene per second key {}\n total multimapper per second key {}", index.info()[0], index.info()[1], index.info()[2], index.info()[3] );

    index.write_index( opts.outpath.to_string() ).unwrap();

    if opts.text{
        index.write_index_txt( opts.outpath.to_string() ).unwrap();
    }

    //index.write_index_txt( opts.outpath.to_string() ).unwrap();
    //eprintln!("THIS IS STILL IN TEST MODE => TEXT INDEX WRITTEN!!! {}",opts.outpath.to_string() );
    

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            eprintln!("finished in {milli} h {min} min {sec} sec {mil} milli sec");
        },
        Err(e) => {println!("Error: {e:?}");}
    }

}