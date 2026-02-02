import pandas as pd
import numpy as np
import sys
import io

def generate_tlens_summary(input_data):
    """
    Analyzes tlens_by_allele data to produce a simplified summary of TL-75 values.
    
    1. Finds the maximum TL-75 (TL_p75) for each unique telomere end (e.g., chr1p, chr1q).
    2. Calculates summary statistics for these maximum TL-75 values,
       both including and excluding sex chromosomes (chrX, chrY).

    Args:
        input_data (str or io.TextIOWrapper): The path to the input TSV file or the file content.
        
    Returns:
        str: A string containing the two-part summary in TSV format.
    """
    
    try:
        # Load the TSV data. pd.read_csv can handle both a file path and a file-like object.
        df = pd.read_csv(input_data, sep='\t')
    except Exception as e:
        # If input_data is a file path and it fails to read, return an error.
        if isinstance(input_data, str):
            # This path is typically taken when running as a script with a filename argument
            # or when a file is explicitly passed.
            sys.stderr.write(f"Error reading file '{input_data}': {e}\n")
            sys.exit(1)
        # This path is taken when reading from stdin or a file-like object
        return f"Error reading data: {e}"

    # Ensure required columns exist
    if '#chr' not in df.columns or 'TL_p75' not in df.columns:
        sys.stderr.write("Error: Input data must contain '#chr' and 'TL_p75' columns.\n")
        sys.exit(1)
    
    # Clean and convert TL_p75 to a numeric type, coercing errors to NaN
    df['TL_p75'] = pd.to_numeric(df['TL_p75'], errors='coerce')
    # Drop any rows where TL_p75 could not be converted
    df.dropna(subset=['TL_p75'], inplace=True)
    
    # --- Part 1: Calculate Maximum TL-75 for each Chromosome End ---

    max_tl_by_end = {}

    # Iterate through the dataframe to handle rows with multiple telomere ends
    for index, row in df.iterrows():
        tl_p75 = int(row['TL_p75'])

        # Split the '#chr' column by comma to get individual telomere ends
        chr_ends = str(row['#chr']).split(',')
        
        for end in chr_ends:
            end = end.strip()
            if not end:
                continue

            # Update the dictionary with the maximum TL_p75 value seen for this end
            max_tl_by_end[end] = max(max_tl_by_end.get(end, -1), tl_p75)

    # Convert the results to a DataFrame for easier processing and sorting
    max_tl_df = pd.DataFrame(max_tl_by_end.items(), columns=['Telomere_End', 'Max_TL_p75'])
    max_tl_df = max_tl_df.sort_values(by='Telomere_End').reset_index(drop=True)
    
    # --- Part 2: Calculate Summary Statistics ---
    
    # All Chromosomes (All unique ends)
    all_tl_values = max_tl_df['Max_TL_p75']
    
    # Autosomes Only (Filter out ends containing 'X' or 'Y', case-insensitive)
    autosomes_df = max_tl_df[
        ~max_tl_df['Telomere_End'].str.contains('X', case=False, na=False) &
        ~max_tl_df['Telomere_End'].str.contains('Y', case=False, na=False)
    ]
    autosomes_tl_values = autosomes_df['Max_TL_p75']

    def calculate_stats(series, name):
        """Calculates min, max, median, mean, count, and stdev."""
        if series.empty:
            return {'Count': 'N/A', 'Min': 'N/A', 'Max': 'N/A', 'Median': 'N/A', 'Mean': 'N/A', 'Stdev_Sample': 'N/A'}
            
        stats = {
            'Count': len(series),
            'Min': series.min(),
            'Max': series.max(),
            'Median': series.median(),
            'Mean': series.mean().round(2),
            # Use ddof=1 for sample standard deviation (stdev)
            'Stdev_Sample': series.std(ddof=1).round(2) 
        }
        
        return stats

    # Get the stats dictionaries
    stats_autosomes = calculate_stats(autosomes_tl_values, 'Autosomes_Only_TL_p75')
    stats_all = calculate_stats(all_tl_values, 'All_Chr_TL_p75')
    
    # Combine the stats into a single DataFrame
    metrics_order = ['Count', 'Min', 'Max', 'Median', 'Mean', 'Stdev_Sample']
    summary_data = {
        'Metric': metrics_order,
        'Autosomes_Only_TL_p75': [stats_autosomes[k] for k in metrics_order],
        'All_Chr_TL_p75': [stats_all[k] for k in metrics_order],
    }
    summary_df = pd.DataFrame(summary_data)
    
    # --- Part 3: Generate Final TSV Output String ---
    
    # Capture Part 1 TSV output
    part1_header = "# Part 1: Maximum TL-75 for each Chromosome End (Max TL-75 per Telomere End)"
    part1_tsv = max_tl_df.to_csv(None, sep='\t', index=False)
    
    # Capture Part 2 TSV output
    part2_header = "\n#\n# Part 2: Overall Summary Statistics on Max TL-75 Values"
    part2_note = "# 'Autosomes_Only' excludes chrX and chrY ends."
    part2_tsv = summary_df.to_csv(None, sep='\t', index=False)

    # Combine all parts
    output = "\n".join([
        "# Simplified TL-75 Summary from tlens_by_allele.tsv",
        "#",
        part1_header,
        "#",
        part1_tsv.strip(),
        part2_header,
        part2_note,
        "#",
        part2_tsv.strip()
    ])
    
    return output

if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        # If not enough arguments, print usage and exit
        sys.stderr.write("Usage: python tel_allele_stats.py <input_tlens_file.tsv> <output_summary_file.tsv>\n")
        sys.exit(1)

    # Get input and output filenames from command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
        
    # Generate the result
    result = generate_tlens_summary(input_file)
    
    # Write the result to the specified output file
    try:
        with open(output_file, 'w') as f:
            f.write(result)
        # Note: In a live environment, you would print a success message, 
        # but in this context, the file is automatically provided.
    except Exception as e:
        sys.stderr.write(f"Error writing to output file '{output_file}': {e}\n")
        sys.exit(1)
