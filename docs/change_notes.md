# Change Notes

## 2025-06-16

- Initial release of the EasySci RNA Pipeline
- Rewrite the `counting_exon_paired_parallel_EasySci.py` to make variables and codes more readable
- Do not require exon gtf file. Normal gtf file is enough.

Current Problems / TODO list:

1. Paired-end reads are counted separately in exon counting, which may cause double counting if R1 and R2 are in the same exon. (Do not affect differential expression analysis)
2. Junction reads are discarded in exon counting, which is a waste of information. If running `samtools view | grep [0-9]N`, 10% of reads are junction reads.
