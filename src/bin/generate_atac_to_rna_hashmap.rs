use std::collections::HashMap;
use std::fs::{File, read_to_string};
use std::io::{self};
use bincode::serialize_into;
use clap::Parser;
use rustody::int_to_str::IntToStr;

#[derive(Parser, Debug)]
#[command(name = "Whitelist Processor")]
#[command(about = "Encodes ATAC and RNA whitelists into a compact binary file.", long_about = None)]
struct Args {
    /// Path to the ATAC whitelist file
    #[arg(short, long)]
    atac: String,

    /// Path to the RNA whitelist file
    #[arg(short, long)]
    rna: String,

    /// Output file path for the binary whitelist map
    #[arg(short, long)]
    output: String,
}

fn process_whitelists(atac_path: &str, rna_path: &str) -> Result<HashMap<u32, u32>, String> {
    let atac_dna = read_to_string(atac_path).map_err(|e| format!("Error reading ATAC file: {}", e))?;
    let rna_dna = read_to_string(rna_path).map_err(|e| format!("Error reading RNA file: {}", e))?;

    let atac_lines: Vec<_> = atac_dna.lines().collect();
    let rna_lines: Vec<_> = rna_dna.lines().collect();

    if atac_lines.len() != rna_lines.len() {
        return Err(format!(
            "ATAC and RNA whitelist files have different lengths: {} vs {}",
            atac_lines.len(),
            rna_lines.len()
        ));
    }

    let mut whitelist_map = HashMap::with_capacity(atac_lines.len());
    let mut int2str = IntToStr::new(b"AAA".to_vec(), 32).unwrap(); // Adjust as needed

    for (atac_line, rna_line) in atac_lines.iter().zip(rna_lines.iter()) {
        let atac_id = int2str.str_to_u32(atac_line);
        let rna_id = int2str.str_to_u32(rna_line);
        whitelist_map.insert(atac_id, rna_id);
    }

    Ok(whitelist_map)
}


fn save_to_binary(map: &HashMap<u32, u32>, output_path: &str) -> io::Result<()> {
    let file = File::create(output_path)?;
    serialize_into(file, map).unwrap();
    Ok(())
}

fn main() {

    let args = Args::parse();

    let combined_map = process_whitelists( &args.atac, &args.rna ).unwrap();

    println!("Saving combined map to binary...");
    save_to_binary(&combined_map, &args.output).expect("Error saving binary file");

    println!("Whitelist binary file created at: {}", &args.output);
}
