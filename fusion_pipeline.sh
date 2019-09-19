cd /home/fuzl/project/amplicon_shengchan_new/
cp 目标序列 ./
for i in  *_1.fq.gz
do
mkdir -p join/${i%%_1.*}_join
o=join/${i%%_1.*}_join
name=${i%%_1*}
fastq-join $i  ${i%%_1.fq.gz}_2.fq.gz -o $o/${i%%_1.fq.gz} >$o/join_${i%%_1.fq.gz}.log && \
mv $o/${i%%_1.fq.gz}join  $o/${name}.fq && \
bwa mem /home/fuzl/project/amplicon_shengchan_new/database_v2/merge_leftright_exon.fa  $o/${name}.fq > $o/${name}.sam  && \
samtools view -Sb $o/${name}.sam |samtools sort >$o/${name}.sorted.bam && \
perl /home/fuzl/script/amplicon_bam_reads_left_right.pl $o/${name}.sorted.bam $o/${name}.sorted.fusion && \
perl /home/fuzl/script/amplicon_left_right_target_filter.pl /home/fuzl/project/amplicon_shengchan_new/database_v2/genename-t.txt \
 $o/${name}.sorted.fusion $o/${name}.sorted.fusion.final &
done
统计
for i in   join/*_join/*.sorted.fusion.final; do echo "$i" ; cut -f 1,2 $i |sort |uniq -c ; done
