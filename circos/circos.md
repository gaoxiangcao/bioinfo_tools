## <center>Circos<center>

### Preparation
- ***extract cds***
  ```bash
  gffread in.gff3 -g ref.fa -x cds.fa
  ```
- ***get pep***
  ```bash
  gffread in.gff3 -g ref.fa -y pep.fa
  ```

### Format Adjustment
- ***title and space***
  ```bash
  cut -d " " -f1 syn.fna | \
  awk '{if($0~">"){print $0}else{printf "%s",$1}}' | \   # $0 refers to the entire current line. $1 refers to the first field (column) of the current line.
  sed 's/>/\n>/g' | grep -A1 -E ">[0-9]" | \
  sed 's/>/>Chr/g' > tair10.genome.fa
  # syn
  sed 's/>NC_000911.1 Synechocystis sp. PCC 6803, complete sequence/>chromsome/g' syn.fna > syn.genome.fa
  
  ```

### chromsome size
- ***chromsome size***
```bash
  genome=syn.genome.fa
  faidx ${genome} -i chromsizes > syn.genome.size
  awk -vFS="\t" -vOFS="\t" '{print "chr","-",$1,"C"NR,"0",$2,"chr"NR}' syn.genome.size > circos.chromosome.txt
  cat circos.chromosome.txt
  ```

### seperation statistic
```bash
  genomeSize=syn.genome.size
  windowSize=10000
  bedtools makewindows -g ${genomeSize} -w ${windowSize} | \
  awk -vFS="\t" -vOFS="\t" '{print $1,$2,$3}' | \
  bedtools sort -i - > ${genomeSize}.10kb
  ```

### gene count 
```bash
  gff3=syn.gff
  genomeSize=syn.genome.size.10kb
  grep '[[:blank:]]gene[[:blank:]]' ${gff3} | \
  awk '{print $1"\t"$4"\t"$5}'| \
  bedtools coverage -a ${genomeSize} -b - | \
  cut -f 1-4 > circos.gene.density.txt
  # color info
  awk -vFS="\t" -vOFS="\t" 'NR==FNR{arr[$3]=$3;brr[$3]=$7}NR!=FNR{if(arr[$1]=$1){print $0,"fill_color="brr[$1]}}' circos.chromosome.txt circos.gene.density.txt > circos.gene.density.color.txt
```

### GC count
```bash
  genome=syn.genome.fa
  genomeSize=syn.genome.size.10kb
  bedtools nuc -fi ${genome} -bed ${genomeSize} | \
  cut -f 1-3,5 | grep -v "#" | \
  awk -vFS="\t" -vOFS="\t" '{print $1,$2,$3,$4*($3-$2)}' > circos.gc.density.txt
  # color info 
  awk -vFS="\t" -vOFS="\t"  'NR==FNR{arr[$3]=$3;brr[$3]=$7}NR!=FNR{if(arr[$1]=$1){print $0,"fill_color="brr[$1]}}' circos.chromosome.txt circos.gc.density.txt > circos.gc.density.color.txt
```

### gene list
```bash
gff3=syn.gff
awk -vFS="\t" -vOFS="\t" '{if($3=="gene"){match($9,/ID=([^;]+)/,a);sub(/ID=/,"",a[0]);print $1,$4,$5,a[0]}}' ${gff3} > ${gff3%%.*}.gene.list.all.txt

grep -f syn_m5c.txt syn.gene.list.all.txt > circos.text.txt
```