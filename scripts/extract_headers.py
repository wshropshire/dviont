def extract_headers_from_fasta(fasta_file):
    """Extract headers (chromosome names) from a multi-FASTA file.
    
    This function assumes that headers are lines that start with '>'.
    It trims everything from the first space to keep only the ID.
    
    Args:
        fasta_file (str): Path to the multi-FASTA file.
    
    Returns:
        list: List of headers (IDs) from the FASTA file, with everything after the first space removed.
    """
    headers = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Remove the '>' symbol and everything after the first space
                headers.append(line[1:].split()[0])  # Keep only up to the first space
    return headers
