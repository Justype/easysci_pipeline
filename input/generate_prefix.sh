#!/bin/bash
N_MIN_RECORDS=5000 # Minimum number of records per sample

if [[ ! -d fastqs ]]; then
  echo "Directory 'fastqs' does not exist. Please create it and place your files inside."
  exit 1
fi

if [[ -f i7_demux_prefix.txt ]]; then
  echo "i7_demux_prefix.txt already exists. No need to run the script again."
  exit 0
fi

# The file names should be in the format: scRNAv3_01_R1.fastq.gz
ls fastqs | cut -d _ -f 1-2 | uniq > i7_demux_prefix.txt

# remove duplicates
sort -u i7_demux_prefix.txt -o i7_demux_prefix.txt.tmp

# Remove IDs with records <= 5000 (5000*4=20000 lines)
rm i7_demux_prefix.txt
for sample in $(cat i7_demux_prefix.txt.tmp); do
  if [[ $(zcat fastqs/${sample}_R1.fastq.gz | grep '^@' | head -n $N_MIN_RECORDS | wc -l) -eq $N_MIN_RECORDS ]]; then
    echo $sample >> i7_demux_prefix.txt.use
  else
    echo "Removing $sample due to insufficient records"
  fi
done

if [[ ! -s i7_demux_prefix.txt.use ]]; then
  echo "No samples with sufficient records found."
  exit 1
fi

# Clean up
mv i7_demux_prefix.txt.use i7_demux_prefix.txt
rm i7_demux_prefix.txt.tmp

# Make sure all files have 4 files _I1, _R1, _R2, _R3
for sample in $(cat i7_demux_prefix.txt); do
  for read in I1 R1 R2 R3; do
    if [[ ! -f fastqs/${sample}_${read}.fastq.gz ]]; then
      echo "Missing file: fastqs/${sample}_${read}.fastq.gz"
    fi
  done
done

