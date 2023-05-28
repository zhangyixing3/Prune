# ALLHIC  Prune in Rust

## install
```
git clone git@github.com:zhangyixing3/Prune.git
cd Prune && cargo build -j 2 --release
```

Finally, the binary file "Prune" can be found in the "Prune/target/release" directory.


## Run
```
$ Prune
new version of prune in Rust

Usage: Prune <COMMAND>

Commands:
  fast  save memory,disk space,and valuable time
  help  Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```
```
$ Prune fast  -h
save memory,disk space,and valuable time

Usage: Prune fast --table <ALLELE> --bam <BAM_F>

Options:
  -i, --table <ALLELE>  Allele.ctg.table
  -b, --bam <BAM_F>     sample.clean.bam
  -h, --help            Print help
  ```


