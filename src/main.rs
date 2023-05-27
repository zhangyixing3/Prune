use rust_htslib::bam::{self,Read, Reader, HeaderView};
use std::collections::{HashMap, HashSet};
mod table_alle;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    author = "Zhang Yixing",
    version = "version 1.1",
    about = "new version of prune in Rust",
    long_about = None
)]
struct Args {
    #[clap(subcommand)]
    command: Subcli,
}
#[derive(Parser, Debug)]
#[allow(non_camel_case_types)]
enum Subcli {
    /// save memory,disk space,and valuable time
    fast{
        /// Allele.ctg.table
        #[arg(short = 'i',long = "table",required=true)]
        allele: String,
        /// sample.clean.bam
        #[arg(short = 'b',long = "bam",required=true)]
        bam:String,
    },
    
}



#[warn(unused_imports)]
#[warn(unused_variables)]
#[warn(unused_imports)]
fn main() {
    let arg = Args::parse();
    let mut alleles: String = String::from("./Allele.ctg.table");
    let mut bamm: String = String::from("./sample.clean.bam");
    match arg.command {
        Subcli::fast { allele, bam } =>{
            alleles = allele;
            bamm = bam;
        },
        _ => eprint!("error!"),
    }
    let alleles = table_alle::read_table(&alleles);
    let mut bamm: Reader = Reader::from_path(&bamm).unwrap();


    // 获取header 信息
    let mut contigs_id: Vec<String> = Vec::new();
    let header = bam::Header::from_template(bamm.header());
    for (_key, records) in header.to_hashmap(){
        for record in records{
            contigs_id.push(record["SN"].clone())
        }
    }
    let mut congtig_index = HashMap::with_capacity(contigs_id.len());
    for (index, value) in contigs_id.iter().enumerate(){
        congtig_index.insert(value, index);
    }

    let mut pairdb:HashMap<u32,HashMap<u32,u32 >> = HashMap::new();
    let mut ctgdb: HashMap<u32, u32> = HashMap::with_capacity(contigs_id.len());
    for read in bamm.rc_records(){
        let read = read.expect("Failure parsing Bam file");
        let mut a :u32 ;
        let mut b :u32 ;
        if read.tid() == read.mtid(){
            continue;
        }else {
            a = read.tid() as u32;
            b = read.mtid() as u32;
            if a > b {
                // Swap values if necessary
                std::mem::swap(&mut a, &mut b);
            }
            let a_entry = pairdb.entry(a).or_insert(HashMap::new());
            let count = a_entry.entry(b).or_insert(0);
            *count += 1;

            ctgdb.entry(a).and_modify(|x |*x += 1).or_insert(1);
            ctgdb.entry(b).and_modify(|x |*x += 1).or_insert(1);


        }
    }
    let index2id: HeaderView = HeaderView::from_header(&header);






    
    let mut removedb: HashMap<u32, HashSet<u32>> = HashMap::new();
    let mut allremovedb: HashMap<u32, HashSet<u32>> = HashMap::new();
    for allele in alleles.clone(){
        for i in 2..allele.len() - 1{
            for j in (i as usize + 1) .. allele.len(){
                // https://rust-bio.github.io/rust-htslib/rust_htslib/bam/record/struct.Record.html#method.tid ref_name to index
                let a_index = index2id.tid(&allele[i].as_ref()).expect("ref name transform error!");
                let b_index = index2id.tid(&allele[j].as_ref()).expect("ref name transform error!");
                // let a_index = &allele[i].name2tid().expect("Invalid reference sequence name");  
                // let b_index = &allele[i].name2tid().expect("Invalid reference sequence name");
                let  a;
                let  b;
                if a_index == b_index{
                    continue;
                }else if a_index < b_index {
                     a = a_index;
                     b = b_index;
                }else {
                    a = b_index;
                    b = a_index;
                }
                removedb.entry(a ).or_insert(HashSet::new()).insert(b );
                allremovedb.entry(a ).or_insert(HashSet::new()).insert(b );
                
            }
        }
    }


    for allele in alleles{   
        for (&key, _value) in &ctgdb{
            let mut retaindb: HashMap<u32, u32> = HashMap::with_capacity(1_usize);

            let mut signal_value = 0;
            for i in 2..allele.len(){
                let ii = index2id.tid(allele[i].as_ref()).expect("ref name transform error!");
                let a: u32 ;
                let b: u32 ;
                if key == ii {
                    continue;
                }else if key > ii  {
                    a = ii;
                    b = key;
                }else {
                    a = key;
                    b = ii;
                }
                if removedb.contains_key(&a) && removedb[&a].contains(&b){
                    continue;
                }
                if !pairdb.contains_key(&a) || !pairdb[&a].contains_key(&b){
                    continue;
                }
                let numr = pairdb[&a][&b];
                match  retaindb.get(&key) {
                    Some(&value) =>{ // 已经有数据，需要比较选出大的值
                        if numr < signal_value {
                            allremovedb.entry(a).or_insert(HashSet::new()).insert(b);
                            // retaindb.entry(i as u32).and_modify(|entry| *entry = b);
                            // signal_value = numr;
                        }else if numr > signal_value {

                            if value > key{
                                allremovedb.entry(key).or_insert(HashSet::new()).insert(value);
                            }else {
                                allremovedb.entry(value).or_insert(HashSet::new()).insert(key);
                            }
                            
                            // allremovedb.entry(a).or_insert(HashSet::new()).insert(b);
                            retaindb.entry(key).and_modify(|entry| *entry = ii);
                            signal_value = numr;
                        }else {
                            print!("{}","equal");
                        }

                    },
                    None => {
                        retaindb.entry(key).and_modify(|entry| *entry = ii);
                        signal_value = numr;
                        
                    },     
                }
            }
        }
    }
    // Iterate through an old BAM file to generate a new pruned BAM file.
    let mut prune_bam = bam::Writer::from_path("./prunning.bam", &header, bam::Format::Bam).expect("creat prunning.bam error !");
    let mut bamm: Reader = Reader::from_path("./merge.bam").unwrap();
    for read in bamm.rc_records(){
        let read = read.expect("Failure parsing Bam file");
        let mut a :u32 ;
        let mut b :u32 ;
        if read.tid() == read.mtid(){
            prune_bam.write(&read).expect("write read error !");
            continue;
        }else {
            a = read.tid() as u32;
            b = read.mtid() as u32;
            if a > b {
                // Swap values if necessary
                std::mem::swap(&mut a, &mut b);
            }
        }
        if allremovedb.get_key_value(&a).is_some() && allremovedb[&a].contains(&b){
            continue;
        }
        prune_bam.write(&read).expect_err("write read error !");
    }
    
}

    




