use clap::Parser;

use regex::Regex;

use this::fast_mapper::FastMapper;
use this::gene::Gene;

use needletail::parse_fastx_file;

//use this::sampleids::SampleIds;
//use this::analysis::

use std::path::PathBuf;
use std::fs;
//use std::path::Path;

use std::time::SystemTime;

use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use flate2::read::GzDecoder;
use std::collections::HashMap;
//use std::collections::HashSet;

use this::mapping_info::MappingInfo;
use this::ofiles::Ofiles;

use std::thread;
use rayon::slice::ParallelSlice;
use rayon::iter::ParallelIterator;

use indicatif::{ProgressStyle, ProgressBar, MultiProgress};

//static EMPTY_VEC: Vec<String> = Vec::new();

// use std::convert::TryInto;

/// Just run a test case for speed reproducability and simplicity

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the gff gene information table
    #[clap(default_value= "testData/testGenes.gtf.gz",short, long)]
    gtf: String,
    /// the fasta genome data
    #[clap( default_value=  "testData/testGenes.fa.gz",short, long)]
    file: String,
    /// the outpath
    #[clap(default_value=  "testData/mapperTest",short, long)]
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




fn process_lines ( lines:&&[String], index: &mut FastMapper ,seq_records: &HashMap< String, Vec<u8>>, chr:&Regex, re_gene_name: &Regex, re_gene_id: &Regex, re_transcript_id: &Regex){

    let mut genes = HashMap::<String, Gene>::new();
    let mut gene_name:String;
    let mut transcript_id:String;
    const COVERED_AREA:usize = 400; // cover 600 bp of the transcript

    for line in lines.iter() {

        let parts: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if parts.len() < 8{
            continue;
        }

        if parts[2] == "transcript"{
            // capture the parts I need using my regexp modules
            if let Some(captures) = re_gene_name.captures( &parts[8].to_string() ){
                gene_name = captures.get(1).unwrap().as_str().to_string();
            }else {
                if let Some(_captures) = re_gene_id.captures( &parts[8].to_string() ){
                    continue; // this likely clutters up the data anyhow.
                }
                else {
                    panic!("I could not identify a gene_name in the attributes {:?}", &parts[8].to_string() );
                }
            }
            // if let Some(captures) = re_gene_id.captures( &parts[8].to_string() ){
            //     gene_id = captures.get(1).unwrap().as_str().to_string();
            // }else {
            //     panic!("I could not identify a gene_id in the attributes {:?}", &parts[8].to_string() );
            // }
            if let Some(captures) = re_transcript_id.captures( &parts[8].to_string() ){
                transcript_id = captures.get(1).unwrap().as_str().to_string();
            }else {
                panic!("I could not identify a transcript_id in the attributes {:?}", &parts[8].to_string() );
            }
            
            // and add a gene
            // pub fn new(chrom:String, start_s:String, end_s:String, sense_strand_s:String, name:String, id:String )
            let gene = Gene::new( parts[0].to_string(),  parts[3].to_string(), parts[4].to_string(), parts[6].to_string(), gene_name.to_string(), vec![transcript_id.to_string()] );
            genes.insert( transcript_id, gene ); 
        }

        if parts[2] == "exon"{
            // capture the parts I need
            if let Some(captures) = re_transcript_id.captures( &parts[8].to_string() ){
                transcript_id = captures.get(1).unwrap().as_str().to_string();
            }else {
                panic!("I could not identify a gene_id in the attributes {:?}", &parts[8].to_string() );
            }
            // and add an exon
            match genes.get_mut( &transcript_id.to_string() ){
                Some(gene) => gene.add_exon( parts[3].to_string(), parts[4].to_string() ),
                None => eprintln!( "ignoring transcript! ({})", transcript_id.to_string())
            }
        }

        if parts[2] == "gene"{
            // here we should start to 'clean out the old ones'!
            let start = match parts[3].parse::<usize>(){
                Ok(v) => v,
                Err(e) => panic!("I could not parse the start of the transcript as usize: {e:?}"),
            };
            let mut to_remove = Vec::new();

            for (k, gene) in &genes {

                if gene.passed(start) {
                    // Do something with the gene, e.g. remove it
                    match seq_records.get( &gene.chrom.to_string() ){
                        Some(seq) => {
                            gene.add_to_index( seq, index, COVERED_AREA );
                            //println!("The genes detected: {:?}", index.names_store );
                        },
                        None => {
                            if chr.is_match ( &gene.chrom.to_string() ){
                                match seq_records.get( &gene.chrom.to_string()[3..] ){
                                    Some(seq) => {
                                        gene.add_to_index( seq, index, COVERED_AREA );
                                        //println!("The genes detected: {:?}", index.names_store );
                                    },
                                    None => {
                                        eprintln!("I do not have the sequence for the chromosome {}", gene.chrom.to_string() );
                                    }
                                }
                            }else {
                                match seq_records.get( &format!("chr{}", &gene.chrom.to_string()) ){
                                    Some(seq) => {
                                        gene.add_to_index( seq, index, COVERED_AREA );
                                        //println!("The genes detected: {:?}", index.names_store );
                                    },
                                    None => {
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
                gene.add_to_index( seq, index, COVERED_AREA );
                //println!("The genes detected: {:?}", index.names_store );
            },
            None => {
                let keys_vec: Vec<_> = seq_records.keys().collect();
                panic!("I do not have the sequence for the chromosome {} {:?}", gene.chrom.to_string(), keys_vec )
            },
        }

    }

}



// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    //// create the report object /////////////////////////////////////
    let opts: Opts = Opts::parse();

    let log_file_str = PathBuf::from(&opts.outpath).join(
        "index_log.txt"
    );

    println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );

    let log_file = match File::create( log_file_str ){
        Ok(file) => file,
        Err(err) => {
            panic!("Error: {err:#?}" );
        }
    };
    let ofile = Ofiles::new( 1, "NOT_USED", "A.gz", "B.gz",  opts.outpath.as_str() );

    let mut report = MappingInfo::new(log_file, 32.0 , 0, ofile );
    report.start_counter();

    //// created the report object /////////////////////////////////////

    let mut kmer_size = opts.gene_kmers;
    if opts.gene_kmers > 32{
        eprintln!("Sorry the max size of the kmers is 32 bp");
        kmer_size = 32;
    }

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

    let mut index = FastMapper::new( kmer_size, 100_000 ); // how many genes do we expect?

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
    
    let gtf = Regex::new(r".*gtf.?g?z?$").unwrap();
    let chr = Regex::new(r"^chr").unwrap();

    let re_gene_name: Regex;
    let re_gene_id: Regex;
    let re_transcript_id : Regex;


    match gtf.is_match( &opts.gtf ){
        true => {
            eprintln!("gtf mode");
            re_gene_name =  Regex::new(r#".* gene_name "([\(\)\w\d\-\._]*)""#).unwrap();
            re_gene_id = Regex::new(r#"gene_id "([\d\w\.]*)";"#).unwrap();
            re_transcript_id = Regex::new(r#"transcript_id "([\d\w\.]*)";"#).unwrap();
        },
        false => {
            eprintln!("gff mode.");
            re_gene_name =  Regex::new(r#".* ?gene_name=([\(\)\w\d\-\._]*); ?"#).unwrap();
            re_gene_id = Regex::new(r#"gene_id=([\d\w\.]*);"#).unwrap();
            re_transcript_id = Regex::new(r#"transcript_id=([\d\w\.]*);"#).unwrap();
        },
    }
    //let mut gene_id:String;
    

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

    let num_threads = num_cpus::get();
    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    pb.set_style(spinner_style);

    let reads_per_chunk = 100_000;
    eprintln!("Starting with data collection");
    let mut good_read_count = 0;
    let max_dim = reads_per_chunk * num_threads;
    let mut lines = Vec::<String>::with_capacity( max_dim );
    let mut batch = 0;

    for line in reader.lines() {
        if good_read_count < max_dim{
            let rec = line.ok().expect("Error reading record.");
            lines.push(rec);
            good_read_count+=1;
        }
        else {
            batch +=1;
            report.stop_file_io_time();
            eprintln!("creating mapper");
            good_read_count = 0;
            let results:Vec<FastMapper> = lines.par_chunks(lines.len() / num_threads + 1) // Split the data into chunks for parallel processing
                .map(|data_split| {
                    // Get the unique identifier for the current thread
                    let _thread_id = thread::current().id();
                
                    // Convert the thread ID to a string for use in filenames or identifiers
                    //let thread_id_str = format!("{:?}",thread_id );
                    //let ofile = Ofiles::new( 1, &("Umapped_with_cellID".to_owned()+&thread_id_str), "R2.fastq.gz", "R1.fastq.gz",  outpath );
                    //let log_file_str = PathBuf::from(outpath).join(
                    //    format!("Mapping_log_{}.txt",thread_id_str )
                    //);
            
                    //let log_file = match File::create( log_file_str ){
                    //        Ok(file) => file,
                    //        Err(err) => {
                    //        panic!("thread {thread_id_str} Error: {err:#?}" );
                    //    }
                    //};
                    let mut idx = FastMapper::new( kmer_size,  reads_per_chunk );
                    // Clone or create a new thread-specific report for each task      
                    let _res = process_lines(&data_split, &mut idx, &seq_records, &chr, &re_gene_name, &re_gene_id, &re_transcript_id );
                    idx

                }) // Analyze each chunk in parallel
            .collect(); // Collect the results into a Vec
            report.stop_multi_processor_time();
            eprintln!("Integrating multicore results");
            for idx in results{
                index.merge(idx);
                //report.merge( &gex.1 );
            }
            report.stop_single_processor_time();
            let (h,m,s,_ms) = MappingInfo::split_duration( report.absolute_start.elapsed().unwrap() );

            eprintln!("For {} sequence regions (x {} steps) we needed {} h {} min and {} sec to process.",max_dim, batch, h, m, s );
            
            //eprintln!("{h} h {m} min {s} sec and {ms} millisec since start");
            eprintln!("{}", report.program_states_string() );

            eprintln!("We created this fast_mapper object:");
            index.eprint();
            eprintln!("Reading up to {} more regions", max_dim);
        }
    }

    eprintln!(" total first keys {}\n total second keys {}\n total single gene per second key {}\n total multimapper per second key {}", index.info()[0], index.info()[1], index.info()[2], index.info()[3] );

    index.write_index( opts.outpath.to_string() ).unwrap();

    if opts.text{
        index.write_index_txt( opts.outpath.to_string() ).unwrap();
    }

    //index.write_index_txt( opts.outpath.to_string() ).unwrap();
    //eprintln!("THIS IS STILL IN TEST MODE => TEXT INDEX WRITTEN!!! {}",opts.outpath.to_string() );
    eprintln!("{}", report.program_states_string() );

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
