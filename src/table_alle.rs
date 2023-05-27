use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;


pub fn read_table(file:&str) -> Vec<Vec<String>>{
    let f = File::open(file).expect("open allele file error");
    let mut lines:Vec<Vec<String>> = Vec::new();
    let reader = BufReader::new(&f);
    for line in reader.lines(){
        let line:Vec<String> = line.unwrap().split('\t').map(|x| x.to_string()).collect();
        if line.len() < 3{
            panic!("The Allele.ctg.table file must have three columns")
        }else if line.len() == 3 {
            continue;
        }else {
            lines.push(line);
        }
    }
    return lines;
}

