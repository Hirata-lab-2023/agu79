#!/bin/bash

# Set output directory
family="agu79"
mkdir -p ${family}/depth ${family}/vcf/HIGH ${family}/vcf/MOD vcf/fin_vcf

# Create necessary directories listed in direc2.txt
for t in `cat direc2.txt`
do
  mkdir -p ${t}
done

# Download FASTQ files from SRA
for t in `cat SRA.txt`
do
  fastq-dump --gzip --split-files --outdir rawdata ${t} &
done
wait

# Trim reads using fastp
for t in `cat rawdata/test.txt`
do
  i=$(grep -n "${t}" rawdata/test.txt | cut -d: -f1)
  fastp --thread 8 \
    --in1 rawdata/${t}/${t}_1.fq.gz \
    --in2 rawdata/${t}/${t}_2.fq.gz \
    --out1 trimed/0${i}_1_trim.fastq.gz \
    --out2 trimed/0${i}_2_trim.fastq.gz \
    --detect_adapter_for_pe &
done
wait

# Align reads to the reference genome using BWA-MEM
for i in {1..6}; do
  bwa mem -t 5 \
    Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa \
    trimed/0${i}_1_trim.fastq.gz trimed/0${i}_2_trim.fastq.gz \
    -o sam/0${i}.sam &
done
wait
rm -f trimed/*

# Convert SAM to BAM, filter out XA tag alignments
for i in {1..6}; do
  samtools view -h sam/0${i}.sam \
    | awk '$17 !~ /XA:/ || $1 ~ /^@/' \
    | samtools view -bS - > bam/0${i}.uniq.bam &
done
wait

# Sort BAM files and add read group information
for i in {1..6}; do
  rm -f sam/0${i}.sam
  samtools sort bam/0${i}.uniq.bam > bam/0${i}.uniq.sort.bam &
done
wait

for i in {1..6}; do
  java -Xmx50G -jar /home/sadamitsu/KS_bin/picard.jar AddOrReplaceReadGroups \
    I=bam/0${i}.uniq.sort.bam \
    O=bam/0${i}.uniq.sort.addG.bam \
    RGID=0${i} RGLB=0${i} RGPL=BGI RGPU=run_barcode RGSM=0${i} &
done
wait

# Index BAM files and mark duplicates
for i in {1..6}; do
  rm -f bam/0${i}.uniq.sort.bam
  samtools index bam/0${i}.uniq.sort.addG.bam &
done
wait

for i in {1..6}; do
  /home/sadamitsu/gatk-4.4.0.0/gatk MarkDuplicates \
    -I bam/0${i}.uniq.sort.addG.bam \
    -O bam/0${i}.uniq.sort.addG_markdup.bam \
    -M bam/0${i}.uniq.sort.addG_markdup_matrics.txt
  samtools index bam/0${i}.uniq.sort.addG_markdup.bam &
done
wait

# Calculate read depth and generate coverage plots
for i in {1..6}; do
  rm -f bam/0${i}.uniq.sort.addG.bam
  samtools coverage bam/0${i}.uniq.sort.addG_markdup.bam > ${family}/depth/0${i}_new_depth.txt &
  tinycov covplot bam/0${i}.uniq.sort.addG_markdup.bam \
    -w 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,MT \
    -s 1000000 -r 1000000 \
    -t ${family}/depth/0${i}_cov.txt -o ${family}/depth/0${i} &
done
wait

# SNP calling, filtering, and annotation
for i in {1..6}; do
  mkdir -p vcf/0${i}
  for chr in `cat chr_mt.txt`; do
    gatk HaplotypeCaller \
      -R Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa \
      -I bam/0${i}.uniq.sort.addG_markdup.bam \
      -L ${chr} \
      -O vcf/0${i}/0${i}_${chr}.vcf &
  done
  wait

  for chr in `cat chr_mt.txt`; do
    bgzip vcf/0${i}/0${i}_${chr}.vcf
    tabix -p vcf vcf/0${i}/0${i}_${chr}.vcf.gz

    gatk SelectVariants \
      -R Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa \
      -V vcf/0${i}/0${i}_${chr}.vcf.gz \
      --select-type-to-include SNP \
      -O vcf/0${i}/0${i}_${chr}.snps.vcf

    gatk VariantFiltration \
      -R Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa \
      -V vcf/0${i}/0${i}_${chr}.snps.vcf \
      -O vcf/0${i}/0${i}_${chr}.snps.MQ60.vcf \
      --filter-expression "QD < 2.0" --filter-name "QD_2.0" \
      --filter-expression "FS > 60.0" --filter-name "FS_60" \
      --filter-expression "MQ < 60.0" --filter-name "MQ_60" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum_-12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_-8.0" \
      --genotype-filter-expression "DP < 4" --genotype-filter-name "lowDP"

    grep -v "^#" vcf/0${i}/0${i}_${chr}.snps.MQ60.vcf | grep "PASS" > vcf/0${i}/0${i}_${chr}_cut.snps.MQ60.vcf
  done

  head -n 100 vcf/0${i}/0${i}_1.snps.MQ60.vcf | grep "^#" > vcf/0${i}/0${i}_1.head.snps.MQ60.vcf
  cat vcf/0${i}/0${i}_1.head.snps.MQ60.vcf vcf/0${i}/0${i}_*_cut.snps.MQ60.vcf > vcf/0${i}/0${i}_all_cat_snps.MQ60.vcf
  bgzip -f vcf/0${i}/0${i}_all_cat_snps.MQ60.vcf
  vcf-sort vcf/0${i}/0${i}_all_cat_snps.MQ60.vcf.gz | bgzip -c > vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz
  tabix -p vcf vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz

  java -jar /home/sadamitsu/snpEff_v4_3t_core/snpEff/snpEff.jar -v \
    -stats vcf/0${i}/0${i}_ano_all.html GRCz11_109 \
    vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz > vcf/0${i}/0${i}_ano_all.vcf

  cat vcf/0${i}/0${i}_ano_all.vcf | java -jar /home/sadamitsu/snpEff_v4_3t_core/snpEff/SnpSift.jar filter "( ANN[*].IMPACT == 'HIGH')" > vcf/0${i}/0${i}_ano_all_HIGH.vcf
  cat vcf/0${i}/0${i}_ano_all.vcf | java -jar /home/sadamitsu/snpEff_v4_3t_core/snpEff/SnpSift.jar filter "( ANN[*].IMPACT == 'MODERATE')" > vcf/0${i}/0${i}_ano_all_MOD.vcf

  bgzip -f vcf/0${i}/0${i}_ano_all.vcf
  tabix -p vcf vcf/0${i}/0${i}_ano_all.vcf.gz
  bgzip -f vcf/0${i}/0${i}_ano_all_HIGH.vcf
  tabix -p vcf vcf/0${i}/0${i}_ano_all_HIGH.vcf.gz
  mv vcf/0${i}/0${i}_ano_all_HIGH.vcf* ${family}/vcf/HIGH/

  bgzip -f vcf/0${i}/0${i}_ano_all_MOD.vcf
  tabix -p vcf vcf/0${i}/0${i}_ano_all_MOD.vcf.gz
  mv vcf/0${i}/0${i}_ano_all_MOD.vcf* ${family}/vcf/MOD/

  cp vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz vcf/fin_vcf/
done
