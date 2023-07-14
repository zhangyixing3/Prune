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
        bam_f:String,
    },
    
}


#[warn(unreachable_patterns)]
#[warn(unused_imports)]
#[warn(unused_variables)]
#[warn(unused_imports)]
fn main() {
    let arg = Args::parse();
    // let mut alleles: String = String::from("./Allele.ctg.table");
    // let mut bam: String = String::from("./sample.clean.bam");
    let (alleles, bam) = match arg.command {
        Subcli::fast { allele, bam_f } =>  (allele, bam_f),
        _ => {
            eprint!("error  command !");
            std::process::exit(1);

        },
    };
    let alleles = table_alle::read_table(&alleles);
    let mut bamm: Reader = Reader::from_path(&bam).unwrap();


    // 获取header 信息, 同时记录bam 文件中ref index 
    let mut contigs_id: Vec<u32> = Vec::new();
    let header = bam::Header::from_template(bamm.header());
    let index2id: HeaderView = HeaderView::from_header(&header);
    for (_key, records) in header.to_hashmap(){
        for record in records{
            if let Some(val) = record.get("SN"){
                let a_index = index2id.tid(val.as_ref()).expect("ref name transform error!");
                contigs_id.push(a_index)
            }
            
        }
    }



    let mut pairdb:HashMap<u32,HashMap<u32,u32 >> = HashMap::new();
    // let mut ctgdb: HashMap<u32, u32> = HashMap::with_capacity(contigs_id.len());
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


        }
    }
    


    

    let mut allremovedb: HashMap<u32, HashSet<u32>> = HashMap::new();
    for allele in alleles.clone(){
        for i in 2..allele.len() - 1{
            for j in (i  + 1) .. allele.len(){
                // https://docs.rs/rust-htslib/latest/rust_htslib/bam/struct.HeaderView.html#method.tid ref_name to index
                // println!("{}-{}",&allele[i],&allele[j]);
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

                allremovedb.entry(a ).or_insert(HashSet::new()).insert(b );
                
            }
        }
    }


    for allele in alleles{   
        for &key in &contigs_id{
            
            let mut retaindb: HashMap<u32, u32> = HashMap::with_capacity(1_usize);

            let mut signal_value = 0;
            for i in 2..allele.len(){
                let ii: u32 = index2id.tid(allele[i].as_ref()).expect("ref name transform error!");
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
                // if allremovedb.contains_key(&a) && allremovedb[&a].contains(&b){
                //     continue;
                // }
                if !pairdb.contains_key(&a) || !pairdb[&a].contains_key(&b){
                    continue;
                }
                let numr = pairdb[&a][&b];
                match  retaindb.get(&key) {
                    Some(&value) =>{ // 已经有数据，需要比较选出大的值
                        if numr <= signal_value {
                            allremovedb.entry(a).or_insert(HashSet::new()).insert(b);
                            // retaindb.entry(i as u32).and_modify(|entry| *entry = b);
                            // signal_value = numr;
                        }else {

                            if value > key{
                                allremovedb.entry(key).or_insert(HashSet::new()).insert(value);
                            }else {
                                allremovedb.entry(value).or_insert(HashSet::new()).insert(key);
                            }
                            
                            // allremovedb.entry(a).or_insert(HashSet::new()).insert(b);
                            retaindb.insert(key, ii);
                            signal_value = numr;
                        }

                    },
                    None => {
                        // retaindb.entry(key).or_insert(ii);
                        retaindb.insert(key, ii);
                        signal_value = numr;
                        
                    },     
                }
                // println!("{:?}",retaindb);
            }
        }
    }
    //     for i in allremovedb.keys(){
    //     println!("{}",i);
    // }

    // Iterate through an old BAM file to generate a new pruned BAM file.
    let mut prune_bam = bam::Writer::from_path("./prunning.bam", &header, bam::Format::Bam).expect("creat prunning.bam error !");
    let mut bamm: Reader = Reader::from_path(&bam).unwrap();
    for read in bamm.rc_records(){
        let read = read.expect("Failure parsing Bam file");
        let mut a :u32 ;
        let mut b :u32 ;
        // println!("{} ------ {}", read.tid(), read.mtid());
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
        prune_bam.write(&read).expect("write read error !");
    }
    
}

    




