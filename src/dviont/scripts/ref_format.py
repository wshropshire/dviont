import os
import logging

def determine_ref_format(ref_file):
    """
    Determine the reference format based on the file extension.
    Args:
        ref_file (str): Path to the reference file.
    Returns:
        str: The reference format ('fasta' or 'genbank').
    """
    # Get the file extension
    file_extension = os.path.splitext(ref_file)[1].lower()  # .fa, .fna, .fasta, .gb, .gbk, .gbbf

    # Check for valid fasta extensions
    if file_extension in ['.fa', '.fna', '.fasta']:
        return 'fasta'

    # Check for valid genbank extensions
    elif file_extension in ['.gb', '.gbk', '.gbff']:
        return 'genbank'

    # If neither, raise an error or set to a default
    else:
        logging.error(f"Unsupported file extension: {file_extension}. Supported extensions are .fa, .fna, .fasta, .gb, .gbk, .gbff.")


    