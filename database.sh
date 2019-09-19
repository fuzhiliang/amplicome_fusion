#重新整数据库
#############

cat /home/zhangsw/reference/annotation/hg19/gencode.v29lift37.annotation.gtf | perl -alne '{next unless $F[2] eq "exon" ;/gene_name \"(.*?)\";.*; exon_number (.*?);/; print "$F[0]\t$F[3]\t$F[4]\t$1\texon$2" }' > total_loc_gene_exon.bed
chr1    11869   12227   DDX11L1 exon1
chr1    12613   12721   DDX11L1 exon2
chr1    13221   14409   DDX11L1 exon3
chr1    12010   12057   DDX11L1 exon1
chr1    12179   12227   DDX11L1 exon2
chr1    12613   12697   DDX11L1 exon3


awk 'NR==FNR{a[$1"\t"$2]}NR>FNR{if ($4"\t"$5 in a){print $0}}' genename-t.txt total_loc_gene_exon.bed >left.bed
awk 'NR==FNR{a[$3"\t"$4]}NR>FNR{if ($4"\t"$5 in a){print $0}}' genename-t.txt total_loc_gene_exon.bed >right.bed

for i in left.bed right.bed
do
bedtools sort -i $i >$i.sorted.bed && bedtools merge -i $i.sorted.bed >$i.merge.bed
done

bedtools intersect -b left.bed.sorted.bed -a left.bed.merge.bed  -wa  -wb|cut -f 1-3,7,8|sort |uniq >left.bed.merge.bed.bed
bedtools intersect -b right.bed.sorted.bed -a right.bed.merge.bed  -wa  -wb|cut -f 1-3,7,8|sort |uniq >right.bed.merge.bed.bed
$ bedtools intersect -b right.bed.sorted.bed -a right.bed.merge.bed  -wa  -wb|cut -f 1-3,7,8|sort |uniq |sed 's/\texon/_exon/'>right.bed.merge.bed.get
$ bedtools intersect -b left.bed.sorted.bed -a left.bed.merge.bed  -wa  -wb|cut -f 1-3,7,8|sort |uniq |sed 's/\texon/_exon/'>left.bed.merge.bed.get
$  bedtools getfasta -fi /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta -bed right.bed.merge.bed  -fo right.bed.merge.fa
$  bedtools getfasta -fi /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta -bed left.bed.merge.bed -fo left.bed.merge.fa

perl  /home/fuzl/script/blast_v2.0/parallel_blast_on_split_query.pl --query right_primer.fa --database right.bed.merge.loc.fa --short  --odir  blast_primer_exon_right
perl  /home/fuzl/script/blast_v2.0/parallel_blast_on_split_query.pl --query left_primer.fa --database left.bed.merge.loc.fa --short  --odir  blast_primer_exon_left

##规范
less right.bed.merge.bed.bed |awk '{print $1"\t"$2"\t"$3"\tright_"$1":"$2"-"$3";"$4"-"$5}'>right.bed.merge.bed.bed.m
less left.bed.merge.bed.bed   |awk '{print $1"\t"$2"\t"$3"\tleft_"$1":"$2"-"$3";"$4"-"$5}'>left.bed.merge.bed.bed.m
bedtools getfasta -fi /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta -bed right.bed.merge.bed.bed.m -name   -fo right.bed.merge.fa
bedtools getfasta -fi /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta -bed left.bed.merge.bed.bed.m -name  -fo left.bed.merge.fa
cd ../database_v2/
cp ../database/left.bed.merge.fa ../database/right.bed.merge.fa ./
cat left.bed.merge.fa right.bed.merge.fa >merge_leftright_exon.fa

