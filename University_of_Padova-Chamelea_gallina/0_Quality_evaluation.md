# Quality evaluation

#### Variables

```bash
Hic_file1="full path hi-c file 1"
Hic_file2="full path hi-c file 2"
HiFi_file="full path HiFi reads file"
Genome_prefix="specify genome prefix"
purge_dups_path="full path to purge_dups bin folder"
```

## Fastqc

```bash
fastqc $Hic_file1 -t 12
fastqc $Hic_file2 -t 12
fastqc $HiFi_file -t 12
```

## Custom python script

```bash
python qstats.py -i $HiFi_file
```
