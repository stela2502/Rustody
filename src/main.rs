use clap::Parser;
//use flate2;
//use flate2::read;
use flate2::Compression;
use flate2::write::GzEncoder;

//use needletail::bitkmer::BitNuclKmer;
//use needletail::{parse_fastx_file, Sequence, FastxReader};
use needletail::{parse_fastx_file, Sequence};
use std::convert::TryInto;

//use serde::Deserialize;
//use std::collections::HashMap;
//use std::io::Write;
//use std::io::BufWriter;
use std::str;
//use anyhow::{Context, Result};

//mod utils;
//use crate::utils::get_all_snps;
//use std::path::Path;
use std::path::PathBuf;
//use std::ffi::OsStr;
use std::fs::File;
use std::fs;

use kmers::naive_impl::Kmer;

use std::collections::BTreeSet;



// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl

/// Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
}

struct Sample {
//    oligo:  String,
    id: usize,
    search: BTreeSet<u64>,
    file1:  GzEncoder<File>, 
    file2:  GzEncoder<File>,
}

impl Sample {
    fn from_description(primer: &[u8], id: usize, sub_len: usize, outpath: &str) -> Self {

        let mut search = BTreeSet::<u64>::new();
        // println!( "split this into {}bp kmers", sub_len );
        // and here comes the fun. Get me all possible sub_len strings in that oligo into the BTreeMap
        for kmer in needletail::kmer::Kmers::new(primer, sub_len.try_into().unwrap() ) {
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            search.insert(km);
        }

        let file1_path = PathBuf::from(outpath).join(format!("{}.fq.gz", id));
        let file2_path = PathBuf::from(outpath).join(format!("{}.fq.gz", id));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let file1 = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
        let file2 = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());
        
        Self {
  //          oligo: String::from_utf8_lossy(primer).into_owned(), // handle errors better
            id,
            search,
            file1,
            file2
        }
    }
}

fn fill_kmer_vec<'a>(seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
   kmer_vec.clear();
   let mut bad = 0;
   for km in seq {
        // I would like to add a try catch here - possibly the '?' works? 
        // if this can not be converted it does not even make sense to keep this info
        for nuc in km{
            if *nuc ==b'N'{
                bad = 1;
            }
        }
        if bad == 0{
            // let s = str::from_utf8(km);
            // println!( "this is the lib: {:?}",  s );
            kmer_vec.push(Kmer::from(km).into_u64());
        }
   } 
}

fn main() -> anyhow::Result<()> {
    // parse the options
    let opts: Opts = Opts::parse();

    // for now just write some tests to the stdout!
    // thanks to Rob I can now skip the stdout part!
    //let stdout = std::io::stdout();
    //let lock = stdout.lock();
    //let buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    let mut samples: Vec<Sample>;// = Vec::with_capacity(12);
    let sub_len = 9;
    // //let File1 = Path::new(outpath).join( Path::new(reads).file_name());
    // let mut File1 = Path::new(&opts.outpath);

    // let mut File2 = File1.join( Path::new(&opts.reads).file_name().unwrap());

    // let file1 = get_writer( &File2 );
    // File2 = File1.join( Path::new(&opts.file).file_name().unwrap() );
    // let file2 = get_writer( &File2);
 
    fs::create_dir_all(&opts.outpath)?;

    if  opts.specie.eq("human") {
        // get all the human sample IDs into this.
        samples = vec![
            Sample::from_description( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", 0, sub_len, &opts.outpath ),    
            Sample::from_description( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", 0, sub_len, &opts.outpath )
        ];

    }
    else if opts.specie.eq("mouse") {
        // and the mouse ones
        samples = vec![
            Sample::from_description( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", 0, sub_len, &opts.outpath ),    
            Sample::from_description( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", 0, sub_len, &opts.outpath ),
            Sample::from_description( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", 0, sub_len, &opts.outpath ),
        ];
    } else {
        println!("Sorry, but I have no primers for species {}", &opts.specie);
        std::process::exit(1)
    }

    let file1_path = PathBuf::from(&opts.outpath).join("ambig.fq.gz");
    let file2_path = PathBuf::from(&opts.outpath).join("ambig.fq.gz");
        
    // need better error handling here too
    //println!( "why does this file break? {}", file1_path.display() );
    let mut file1_ambig_out = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
    let mut file2_ambig_out = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());

    // for now, we're assuming FASTQ and not FASTA.
    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    let mut kmer_vec = Vec::<u64>::with_capacity(60);
    let mut unknown = 0;

    while let Some(record2) = readefile.next() {
        if let Some(record1) = readereads.next() {
            let seqrec = record2.expect("invalid record");
            let seqrec1 = record1.expect("invalid record");
            //let seq = seqrec.seq().into_owned();

            let norm_seq = seqrec.normalize(true);

            // create 9mers of the R2 read and check which of the 12 sample ids matches best:
            // println!( "split R2 into {}bp kmers", sub_len );

            let kmers = norm_seq.kmers(sub_len.try_into().unwrap());
            //let kmers = seqrec.kmers(sub_len.try_into().unwrap());
            fill_kmer_vec(kmers, &mut kmer_vec);

            let mut res = vec![0; 12]; //Vec::with_capacity(12);
            let mut max_value = 0;

            for i in 0..samples.len(){
                for s in &kmer_vec {
                    if samples[i].search.contains(&s){
                        //println!( "kmer matches to sample {}", samples[i].id );
                        res[i] +=1;
                    }
                }
                if res[i] > max_value{
                    max_value = res[i];
                }
            }

            let mut z = 0;
            let mut id = 0;
            if max_value > 2 {
                for i in 0..res.len(){
                    if res[i] == max_value {
                        id = i;
                        z += 1;
                    }
                }
            }
            //println!( "we have a match to sample {} with a max value of {} and {} samples reaching this value ",id, max_value,z );
            if z == 1 {
                samples[id].id += 1;
                seqrec1.write(&mut samples[id].file1, None)?;
                seqrec.write(&mut samples[id].file2, None)?;
            } else {
                unknown += 1;
                seqrec1.write(&mut file1_ambig_out, None)?;
                seqrec.write(&mut file2_ambig_out, None)?;
            }

        } else {
            anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
        }
    }
    
    println!( "collected sample info:");
    for i in 0..samples.len(){
        println!( "    sample {}: {} reads", i, samples[i].id );
        samples[i].file1.try_finish()?;
        samples[i].file2.try_finish()?;
    }
    println!( "      unknown: {} reads", unknown );
    file1_ambig_out.try_finish()?;
    file2_ambig_out.try_finish()?;
    


    Ok(())
}
