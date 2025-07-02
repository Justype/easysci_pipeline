#!/bin/bash

# Backup exon counting results by given name

if [ -z "$1" ]; then
  echo "Usage: $0 <backup_name>"
  exit 1
fi

backup_name="$1"
final_dir="output/final"
human_backup_dir="$final_dir/human_exon_$backup_name"
mouse_backup_dir="$final_dir/mouse_exon_$backup_name"

printf "Rename existing $final_dir/exons and remove output/{human,mouse}/exon. Are you sure? [y/N]: "
read -r confirmation
if [[ "$confirmation" != "y" && "$confirmation" != "Y" ]]; then
    echo "Renaming cancelled."
    exit 0
fi

if [ -d "$human_backup_dir" ]; then
  echo "Human exon backup directory already exists: $human_backup_dir"
  exit 1
else
  if [ -d "$final_dir/human_exon" ]; then
    echo "Renaming existing human_exon directory to $human_backup_dir"
    mv "$final_dir/human_exon" "$human_backup_dir"
  else
    echo "No existing human_exon directory to rename."
  fi
fi

if [ -d "$mouse_backup_dir" ]; then
  echo "Mouse exon backup directory already exists: $mouse_backup_dir"
  exit 1
else
  if [ -d "$final_dir/mouse_exon" ]; then
    echo "Renaming existing mouse_exon directory to $mouse_backup_dir"
    mv "$final_dir/mouse_exon" "$mouse_backup_dir"
  else
    echo "No existing mouse_exon directory to rename."
  fi
fi

# Remove the current exon directories
if [ -d "output/human/exon" ]; then
  echo "Removing output/human/exon"
  rm -r output/human/exon
fi
if [ -d "output/mouse/exon" ]; then
  echo "Removing output/mouse/exon"
  rm -r output/mouse/exon
fi
