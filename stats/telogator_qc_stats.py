import pandas as pd
import argparse
import re
import sys
import os

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def main():
    parser = argparse.ArgumentParser(description="Generate a chromosome ends report from telomere length data.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input tlens_by_allele.tsv file")
    parser.add_argument("-f", "--fai", required=True, help="Path to the genome .fai file to define expected chromosomes")
    parser.add_argument("-o", "--output", default="telomere_report.txt", help="Path for the output report file (default: telomere_report.txt)")
    
    # Show help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    
    # Check for file existence
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
    if not os.path.exists(args.fai):
        print(f"Error: FAI file '{args.fai}' not found.")
        sys.exit(1)

    # 1. Process FAI to determine all expected chromosome ends (p and q)
    expected_ends = []
    try:
        with open(args.fai, 'r') as f:
            for line in f:
                if line.strip():
                    chrom = line.split('\t')[0]
                    expected_ends.append(f"{chrom}p")
                    expected_ends.append(f"{chrom}q")
    except Exception as e:
        print(f"Error reading FAI file: {e}")
        sys.exit(1)
        
    # 2. Process TSV for found/ambiguous ends
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
        
    found_ends = set()
    ambiguous_set = set()
    total_count = 0
    unambiguous_count = 0
    
    # Ensure correct column name handling
    chr_col = '#chr' if '#chr' in df.columns else 'chr'
    if chr_col not in df.columns:
        print(f"Error: Column '#chr' (or 'chr') not found in {args.input}")
        sys.exit(1)
    
    for entry in df[chr_col]:
        # Handle comma-separated multi-mappings
        parts = [p.strip() for p in str(entry).split(',')]
        total_count += len(parts)
        
        if len(parts) > 1:
            for p in parts:
                ambiguous_set.add(p)
        else:
            unambiguous_count += 1
            
        for p in parts:
            found_ends.add(p)
            
    # 3. Calculate missing ends (Expected - Found)
    missing_ends = sorted([e for e in expected_ends if e not in found_ends], key=natural_sort_key)
    ambiguous_ends_sorted = sorted(list(ambiguous_set), key=natural_sort_key)
    
    # 4. Format the output table
    output_lines = []
    output_lines.append("Chromosome ends")
    output_lines.append("")
    output_lines.append(f"{'unambiguous':<15} {'total':<10}")
    output_lines.append(f"{unambiguous_count:<15} {total_count:<10}")
    output_lines.append("")
    output_lines.append(f"{'ambiguous ends':<23} {'missing end':<20}")
    
    max_rows = max(len(ambiguous_ends_sorted), len(missing_ends))
    
    if max_rows == 0:
        output_lines.append(f"{'none':<23} {'none':<20}")
    else:
        for i in range(max_rows):
            ambig = ambiguous_ends_sorted[i] if i < len(ambiguous_ends_sorted) else ""
            if not missing_ends:
                miss = "none" if i == 0 else ""
            else:
                miss = missing_ends[i] if i < len(missing_ends) else ""
            output_lines.append(f"{ambig:<23} {miss:<20}")
            
    # 5. Write to file
    try:
        with open(args.output, 'w') as out:
            out.write("\n".join(output_lines) + "\n")
        print(f"Report successfully saved to: {args.output}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
