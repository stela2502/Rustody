use clap::Parser;
use std::env;

//use regex::Regex;

use rustody::genes_mapper::GenesMapper;
use rustody::mapping_info::MappingInfo;

use std::time::SystemTime;
use std::fs;


#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the fasta file (should be gzipped!)
    #[clap(short, long)]
    file: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// how many threads to use to analyze this (default 1)
    #[clap( short, long)]
    num_threads: Option<usize>,
    /// How many genes may each 16bp fragment link to (default 10)
    #[clap( short, long)]
    max_links: Option<usize>,
}

// Function to check if a file exists
fn check_file_existence(file_path: &str, option: &str, errors: &mut Vec<String>) {
    if fs::metadata(file_path).is_err() {
        errors.push(format!("Option {option} - File not found: {file_path}"));
    }
}


// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    //// create the report object /////////////////////////////////////
    let opts: Opts = Opts::parse();

    // Check if each file exists
    let mut errors = Vec::new();

    check_file_existence(&opts.file, "file", &mut errors);

    let num_threads = match &opts.num_threads{
        Some(n) => *n,
        None => 1,
    };
    let max_links = match &opts.max_links{
        Some(n) => *n,
        None => 10,
    };

    // If there are errors, print them and exit
    if !errors.is_empty() {
        eprintln!("Error: Some files do not exist:");
        for error in &errors {
            eprintln!("{}", error);
        }
        std::process::exit(1);
    }

    if fs::metadata(&opts.outpath).is_ok() {
        if let Err(err) = fs::remove_dir_all(&opts.outpath) {
            eprintln!("Error old index directory: {}", err);
            std::process::exit(1);
        } else {
            println!("Old index directory removed successfully!");
        }

        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory: {}", err);
            std::process::exit(1);
        } else {
            println!("New index directory created successfully!");
        }
    }else if let Err(err) = fs::create_dir_all(&opts.outpath) {
        eprintln!("Error creating directory: {}", err);
        std::process::exit(1);
    } else {
        println!("New index directory created successfully!");
    }

    /*let log_file_str = PathBuf::from(&opts.outpath).join(
        "index_log.txt"
    );*/


    let mut report = MappingInfo::new( None, 32.0 , 0, None);
    report.start_counter();

    //// created the report object /////////////////////////////////////

    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

    let mut genes = GenesMapper::new(0);

    match genes.from_fastq( Some(opts.file) ){
        Ok(()) => {
            // all is fine
        },
        Err(e) => {
            eprintln!("{e:?}")
        }
    }

    report.stop_file_io_time();

    genes.filter( max_links, num_threads, );

    report.stop_multi_processor_time();

    let _ = genes.write_index( &opts.outpath );

    report.stop_file_io_time();

    println!("{}", report.program_states_string() );

    let prog_name = env::current_exe()
        .ok()
        .and_then(|p| p.file_name()?.to_str().map(String::from)) // Convert to `String`
        .unwrap_or_else(|| "unknown program".into()); // Provide default

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("{prog_name} finished in {milli}h {min}min {sec} sec {mil}milli sec\n" );},
       Err(e) => {println!("Error: {e:?}");}
    }


}