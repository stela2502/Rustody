/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use std::collections::HashSet;

use kmers::naive_impl::Kmer;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;

// #[derive(Debug)]
// pub struct Info {
// //    pub cells: CellIds10x,
//     pub id: u64,
//     pub name: std::string::String
// }


// impl  Info{
//         pub fn new( id: u64, name: std::string::String )-> Self {
//             let loc_name = name.clone();
//             Self {
//                 id,
//                 name: loc_name,
//             }
//         }
// }

/// GeneIds harbors the antibody tags
/// these sequences have (to my knowmledge) 15 bp os length
/// but that can be different, too. 
/// Hence we store here 
/// kmers       : the search object
/// seq_len     : the length of the sequences (10x oversequences them)
/// kmer_size   : the length of the kmers
/// names       : a hashset for the gene names
/// bad_entries : a hash to save bad entries (repetetive ones)
pub struct GeneIds{    
    pub kmers: BTreeMap<u64, usize>, // the search map with kamer u64 reps.
    pub seq_len: usize, // size of the sequence that has been split into kmers
    kmer_size: usize, // size of the kmers
    names_store: BTreeMap< usize, std::string::String>,
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    bad_entries: HashSet<u64>, // non unique u64 values that will not be recoreded.
    pub good_entries: HashSet<u64>, // upon export as sparse matrix it has to be checked if a gene has a value
    pub max_id: usize,// hope I get the ids right this way...
    unchecked: bool, //if add_unchecked was used results should always be as get_unchecked
}

// here the functions
impl GeneIds{
    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new(kmer_size: usize )-> Self {
        let kmers = BTreeMap::<u64, usize>::new();
        let names = BTreeMap::<std::string::String, usize>::new();
        let names_store = BTreeMap::< usize, std::string::String>::new();
        let names4sparse = BTreeMap::<std::string::String, usize>::new();
        let bad_entries = HashSet::<u64>::new();
        let good_entries = HashSet::<u64>::new();
        let seq_len = 0;
        let max_id = 0;
        let unchecked = false;
        Self {
            kmers,
            seq_len,
            kmer_size: kmer_size,
            names_store,
            names,
            names4sparse,
            bad_entries,
            good_entries,
            max_id,
            unchecked
        }
    }

    pub fn add(&mut self, seq: &[u8], name: std::string::String ){
        
        if seq.len() > self.seq_len{
            self.seq_len = seq.len() 
        }

        let mut checker = BTreeMap::<u8, usize>::new();
        let mut total = 0;

        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            checker.clear();
            for nuc in kmer {
                match checker.get_mut( nuc ){
                    Some( map ) => *map += 1,
                    None => {
                        match checker.insert( nuc.clone(), 1){
                            Some(_) => (),
                            None => (),
                        }
                    }
                }
                if *nuc ==b'N'{
                    continue;
                }
            }
            if checker.len() < 3{
                //println!( "kmer for gene {} is too simple/not enough diff nucs): {:?}", name, std::str::from_utf8(kmer)  );
                continue;
            }
            for ( _key, value ) in checker.iter(){
                if *value as f32 / self.kmer_size as f32 > 0.6 {
                    //println!( "kmer for gene {} is too simple/too many nucs same: {:?}", name, std::str::from_utf8(kmer)  );
                    continue;
                } 
            }

            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            if self.bad_entries.contains( &km ){
                continue
            }
            if self.kmers.contains_key ( &km ){
                self.bad_entries.insert( km.clone() );
                self.kmers.remove( &km );
            }else {
                //let info = Info::new(km, name.clone() );
                if ! self.names.contains_key( &name ){
                    self.names.insert( name.clone(), self.max_id );
                    self.names_store.insert( self.max_id, name.clone() );
                    self.max_id += 1;
                }
                //println!("I insert a kmer for id {}", self.max_id-1 );
                self.kmers.insert(km, self.max_id-1 );
                total +=1;
            }
        }
        if total < 3{
            eprintln!( "Sequence for gene {} has too little OK kmers to match!", name);
        }
        //println!( "{} kmers for gene {}", total, name );

    }


    // TTTATTCATACAGTGCATCGAAATTGGTGTTGCGGTATAGGGACTCCAAAAAAAAAAAAAAAAAAAAAAAT
    //TAAAACACTTCTGTGCATCGAAATTGGTGTTGACGGTATAGGGACTCCAAAAAAAAAAAAAAAAAAAAAAA

    pub fn add_unchecked(&mut self, seq: &[u8], name: std::string::String ){
        
        if seq.len() > self.seq_len{
            self.seq_len = seq.len() 
        }
        self.unchecked = true;
        let mut checker = BTreeMap::<u8, usize>::new();
        let mut total = 0;

        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            checker.clear();
            for nuc in kmer {
                match checker.get_mut( nuc ){
                    Some( map ) => *map += 1,
                    None => {
                        match checker.insert( nuc.clone(), 1){
                            Some(_) => (),
                            None => (),
                        }
                    }
                }
                if *nuc ==b'N'{
                    continue;
                }
            }
            if checker.len() < 3{
                //println!( "kmer for gene {} is too simple/not enough diff nucs): {:?}", name, std::str::from_utf8(kmer)  );
                continue;
            }
            for ( _key, value ) in checker.iter(){
                if *value as f32 / self.kmer_size as f32 > 0.6 {
                    //println!( "kmer for gene {} is too simple/too many nucs same: {:?}", name, std::str::from_utf8(kmer)  );
                    continue;
                } 
            }

            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            if self.bad_entries.contains( &km ){
                continue
            }
            if self.kmers.contains_key ( &km ){
                self.bad_entries.insert( km.clone() );
                self.kmers.remove( &km );
            }else {
                //let info = Info::new(km, name.clone() );
                if ! self.names.contains_key( &name ){
                    self.names.insert( name.clone(), self.max_id );
                    self.names_store.insert( self.max_id, name.clone() );
                    self.max_id += 1;
                }
                //println!("I insert a kmer for id {}", self.max_id-1 );
                self.kmers.insert(km, self.max_id-1 );
                total +=1;
            }
            //println!("I insert a kmer for id {}", self.max_id-1 );
            self.kmers.insert(km, self.max_id-1 );
            total +=1;
        }
        if total < 3{
            eprintln!( "Sequence for gene {} has too little OK kmers to match!", name);
        }
        //println!( "{} kmers for gene {}", total, name );

    }
    pub fn get_id( &mut self, name: String ) -> usize{
        let id = match self.names.get( &name ) {
            Some( id ) => id,
            None => panic!("Gene {} not defined in the GeneID object", name),
        };
        return *id ;
    }

    pub fn get_name( &self, id:usize) -> String{
        let name = match self.names_store.get( &id ){
            Some( na ) => na,
            None => panic!("GeneID {} not defined in the GeneID object", id),
        };
        return name.to_string()
    }

    pub fn get(&mut self, seq: &[u8] ) -> Option< usize >{
        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        if self.unchecked{
            return self.get_unchecked( seq );
        }
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        //let report = self.get_id("2810417H13Rik".to_string());

        let mut ret:Option<usize> = None;

        if kmer_vec.len() == 0 {
            //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
            return ret
        }  
        let mut sums = vec![0 ;self.names.len()];
        let mut max = 0;

        for km in kmer_vec{
            //println!( "searching for kmer {}", km);
            match self.kmers.get(&km){
                Some(c1) => {
                    //println!("And got a match: {}", c1);
                    sums[*c1] += 1;
                    if max < sums[*c1]{
                        //println!("the new max is {}", max);
                        max =  sums[*c1];
                        if max > 4{
                            // 2 hits gave me really really strange results.
                            // need to check for 3 again.
                            break // 2 unique hits should be enough
                        }
                    };
                }
                None => ()
            };
        }
        if max >4 {
            for i in 0..sums.len(){
                if sums[i] == max{
                    //println!("Now the ret hould have the value {} resp {:?}", i, Some(i));
                    if ! ret.is_none() {
                        //eprintln!("I have two genes matching with max value of {}: {:?} and {}", max, ret, i);
                        ret =None;
                        return ret
                    }
                    ret = Some(i);
                }
            }
        }
        //println!("return geneid {:?}", ret);
        return ret
    }

    pub fn get_unchecked(&mut self, seq: &[u8] ) -> Option< usize >{
        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut ret:Option<usize> = None;

        if kmer_vec.len() == 0 {
            //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
            return ret
        }  
        let mut sums = vec![0 ;self.names.len()];
        let mut max = 0;

        for km in kmer_vec{
            //println!( "searching for kmer {}", km);
            match self.kmers.get(&km){
                Some(c1) => {
                    //println!("And got a match: {}", c1);
                    sums[*c1] += 1;
                    if max < sums[*c1]{
                        //println!("the new max is {}", max);
                        max =  sums[*c1];
                    };
                }
                None => ()
            };
        }

        if max >3 {
            for i in 0..sums.len(){
                if sums[i] == max{
                    if ! ret.is_none() {
                        //eprintln!("I have two genes matching with max value of {}: {:?} and {}", max, ret, i);
                        ret =None;
                        return ret
                    }
                    //println!("max = {} -> now the ret hould have the value {:?}", max, Some(i));
                    ret = Some(i);
                    //break;
                }
            }
        }
        // else {
        //     eprintln!("The max was {} => no (good) gene found", max);
        // }
        //println!("return geneid {:?}", ret);
        return ret
    }
    // pub fn to_ids( &self,  ret:&mut Vec<Info> )  {
    //     ret.clear();
    //     for (_i, obj) in &self.kmers {
    //         ret.push(*obj );
    //     }
    // }

    pub fn to_header( &self ) -> std::string::String {
        let mut ret= Vec::<std::string::String>::with_capacity( self.names.len() +2 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for (obj, _id) in &self.names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( format!( "{}", obj)) ;
        }
        ret.push("Most likely name".to_string());
        ret.push("Faction total".to_string());
        return "CellID\t".to_owned()+&ret.join("\t")
    }

    pub fn to_header_n( &self, names: &Vec<String> ) -> std::string::String {
        let mut ret= Vec::<std::string::String>::with_capacity( names.len() +2 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for name in names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( format!( "{}", name)) ;
        }
        ret.push("AsignedSampleName".to_string());
        ret.push("FractionTotal".to_string());
        return "CellID\t".to_owned()+&ret.join("\t")
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


#[cfg(test)]
mod tests {

    use crate::geneids::GeneIds;
     #[test]
    fn check_geneids() {
        let mut genes = GeneIds::new( 7 );

        let mut geneid = 0;
        
        genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
        genes.names4sparse.insert( "Gene1".to_string(), geneid );
        genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
        genes.names4sparse.insert( "Gene1".to_string(), geneid );
        // the next two should not be in the output
        genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
        genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

        geneid = genes.get_id( "Gene1".to_string() );
        assert_eq!( geneid,  0 ); 

        geneid = genes.get_id( "Gene3".to_string() );
        assert_eq!( geneid,  2 ); 
            
    }
}