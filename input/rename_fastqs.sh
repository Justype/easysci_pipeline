#!/bin/bash

################################################################
# Script to rename fastq files and create i7_demux_prefix.txt
# - Example: scRNAv3_P7-1_01_S01_R1_001.fastq.gz
# - Renamed to: scRNAv3_01_R1.fastq.gz
#
# If i7_demux_prefix.txt already exists, the script will exit.
#
# Make sure you have the 'orad' command available to convert fastq.ora files.
################################################################


# SETTINGS
N_MIN_RECORDS=5000 # Minimum number of records per sample

if [[ ! -d fastqs ]]; then
  echo "Directory 'fastqs' does not exist. Please create it and place your files inside."
  exit 1
fi

if [[ -f i7_demux_prefix.txt ]]; then
  echo "i7_demux_prefix.txt already exists. No need to run the script again."
  exit 0
fi

for file in fastqs/*fastq.ora; do
  orad "$file" # orad will convert fastq.ora to fastq.gz
done

for file in fastqs/*fastq.gz; do
  if [[ $file =~ .*?_[RI][1234]_001.fastq.gz ]]; then
    printf "Processing $file\n"
    # Example name: scRNAv3_P7-1_01_S01_R1_001.fastq.gz
    # We want to rename it to: scRNAv3_01_R1.fastq.gz
    new_name=$(basename $file | cut -d _ -f 1,3,5).fastq.gz
    mv $file fastqs/$new_name

    echo $new_name | cut -d _ -f 1,2 >> i7_demux_prefix.txt

    # also edit .md5
    if [[ ! -f $file.md5 ]]; then
      echo "No md5 file found for $file"
      continue
    fi
    mv $file.md5 fastqs/$new_name.md5
    sed -i "s/$(basename $file)/$new_name/" fastqs/$new_name.md5
  fi
done

if [[ ! -f i7_demux_prefix.txt ]]; then
  echo "i7_demux_prefix.txt not found. Exiting."
  echo "Please ensure that fastqs directory contains valid fastq files."
  exit 1
fi

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

