# Change Notes

<a target="_blank" href="https://colab.research.google.com/github/Justype/easysci_pipeline/blob/main/docs/ipynbs/EasySci_issues.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a> EasySci Exon Counting Issues

## 2025-06-18

- Add different exon counting methods:
  - `balanced`: treats R1 and R2 as a group, and counts reads in the same exon only once.
  - `junction`: uses both exon and junction reads and only counts once if both blocks are in the same exon (blocks are separated by CIGAR operator N).
- Change `mem_mb` to `mem_mib` in `config.yaml`.

## 2025-06-16

- Initial release of the EasySci RNA Pipeline
- Rewrite the `counting_exon_paired_parallel_EasySci.py` to make variables and codes more readable
- Do not require exon gtf file. Normal gtf file is enough.

Current Problems / TODO list:

1. Paired-end reads are counted separately in exon counting, which may cause double counting if R1 and R2 are in the same exon. (Do not affect differential expression analysis)
2. Junction reads are discarded in exon counting, which is a waste of information. If running `samtools view | grep [0-9]N`, 5-15% of reads are junction reads.
