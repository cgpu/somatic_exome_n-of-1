#! /usr/bin/env bash

echo "##fileformat=VCFv4.1" > "candidates.vcf"
samtools view -H $2 | grep @SQ | while read LINE; do
  CHR=$(echo $LINE | perl -ane '$F[1]=~s/SN://; print $F[1];')
  LEN=$(echo $LINE | perl -ane '$F[2]=~s/LN://; print $F[2];')
  echo "##contig=<ID=${CHR},length=${LEN}>"
done >> "candidates.vcf"
echo -e "#CHROM\tPOS\tREF\tALT" >> "candidates.vcf"
for VCF in *.vcf; do
  grep -v '#' $VCF | cut -f 1,2,4,5;
done | sort -V | uniq >> "candidates.vcf"
