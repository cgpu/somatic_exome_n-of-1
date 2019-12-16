# Exome N-of-1 Analysis
## Call SNV and CNA from multiple tumour sample exome data
### How to Run
#### Reference Generation
The NextFlow pipelines require references (fasta, indexes, etc.) to run. With dependencies NextFlow and Singularity installed, these can be generated using our DNAseq_references Nextflow
```
git clone https://github.com/brucemoran/DNAseq_references
cd DNAseq_references
nextflow run DNAseq_references.GRCh37.simg.nf \
  --exomebed "/local/full/path/to/exome.bed" \
  --cosmicCGCtsv "/local/full/path/to/COSMIC_Census.tsv" \
  -c DNAseq_references.GRCh37.simg.nextflow.config
```

#### Main pipeline
##### Exome N-of-1
To then run exome_n-of-1 pipeline:
```
nextflow run main.nf \
  --sampleCsv "sample.csv" \
  --refDir "GRCh37" \
  --includeOrder "tumour_A,tumour_B" \
  --germline
```
sample.csv has format:
```
type,sampleID,meta,read1,read2
germline,germ1,"germline sample",/full/path/to/germ1.R1.fastq.gz,/full/path/to/germ1.R2.fastq.gz
somatic,soma1,"somatic sample 1",/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
somatic,soma2,"somatic sample 2",/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
```

'meta' is used for reporting output, where 'sampleID' may include clinically sensitive patient information

N.B that headers of csv files must match exactly, and you must have only one germline/normal sample; please give type as 'germline' for the sake of it being correctly identified.
