
def convert_fai_to_indexes(fai_file):

    import sys

    # Initialize the variables
    # The desired structure for LEXICO_2_IND maps chromosome name to a numeric index.
    lexico_2_ind = {}
    # T2T_CHROMSIZE maps chromosome name to its length (size).
    t2t_chromsize = {}

    chrom_digits = []

    # Keep track of the index for LEXICO_2_IND
    # This index will be assigned sequentially based on the order in the FAI file.
    current_index = 1

    try:
        with open(fai_file, 'r') as f:
            # Read the file line by line
            for line in f:
                # An FAI file is tab-separated. Split the line into fields.
                fields = line.strip().split('\t')

                if len(fields) >= 2:
                    chrom_name = fields[0]
                    chr_digit = chrom_name.replace("chr","")
                    chrom_digits.append(str(chr_digit))
                    try:
                        chrom_size = int(fields[1])
                    except ValueError:
                        # Skip lines where the size is not a valid integer
                        print(f"Warning: Could not convert size for '{chrom_name}'. Skipping.", file=sys.stderr)
                        continue

                    # Populate T2T_CHROMSIZE
                    t2t_chromsize[chrom_name] = chrom_size

                    # Populate LEXICO_2_IND with the sequential index
                    lexico_2_ind[chrom_name] = current_index

                    # Increment the index for the next chromosome
                    current_index += 1
    except FileNotFoundError:
        print(f"Error: FAI file not found at '{fai_file}'. Please check the file path.", file=sys.stderr)
        # Exit the script or handle the error as appropriate
        sys.exit(1)
    if "chrU" not in lexico_2_ind:
        lexico_2_ind["chrU"] = current_index
    return t2t_chromsize, lexico_2_ind, chrom_digits
