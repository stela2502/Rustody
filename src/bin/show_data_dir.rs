
//used for the install process
fn main() {
    let path = bam_tide::get_data_dir();
    println!("{}", path.display());
}