#!/usr/bin/env python3

import pandas as pd
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge count logs for human and mouse samples.")
    parser.add_argument("-i", "--input", type=str, default="output/logs", help="Input directory containing log files.")
    parser.add_argument("-o", "--output", type=str, default="output/final", help="Output directory for merged stats.")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Input directory {args.input} does not exist. Please check the path.")
        return

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    human_input_path = os.path.join(args.input, "human")
    mouse_input_path = os.path.join(args.input, "mouse")

    human_output_path = os.path.join(args.output, "human_stats.csv")
    mouse_output_path = os.path.join(args.output, "mouse_stats.csv")

    # Check the input logs and get the prefix by the name {prefix}.log

    if os.path.exists(human_input_path):
        human_stats = pd.DataFrame()

        for human_file in sorted(os.listdir(human_input_path)):
            if human_file.endswith(".log"):
                i7_prefix = human_file[:-4]  # Remove the .log extension
                stats = pd.read_csv(os.path.join(human_input_path, human_file), sep="\t", index_col=1, header=None)
                stats.columns = [i7_prefix]
                human_stats = pd.concat([human_stats, stats], axis=1)  # column-wise concatenation
        
        human_stats = human_stats.fillna(-1).astype(int).transpose()
        human_stats.to_csv(human_output_path)
    
    if os.path.exists(mouse_input_path):
        mouse_stats = pd.DataFrame()

        for mouse_file in sorted(os.listdir(mouse_input_path)):
            if mouse_file.endswith(".log"):
                i7_prefix = mouse_file[:-4]
                stats = pd.read_csv(os.path.join(mouse_input_path, mouse_file), sep="\t", index_col=1, header=None)
                stats.columns = [i7_prefix]
                mouse_stats = pd.concat([mouse_stats, stats], axis=1)
        mouse_stats = mouse_stats.fillna(-1).astype(int).transpose()
        mouse_stats.to_csv(mouse_output_path)

if __name__ == "__main__":
    main()