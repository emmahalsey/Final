import sys

def validate_sequence(sequence, k):
    """
    Checks whether a DNA sequence is valid for k-mer analysis.

    A valid sequence:
    - Must be at least length k
    - Must contain only valid DNA characters (A, C, G, T)

    Args:
        sequence (str): The DNA sequence to validate
        k (int): Length of k-mers

    Returns:
        bool: True if valid, False otherwise
    """

    # Sequence must be at least length k
    if len(sequence) < k:
        return False

    # Check that each character is a valid DNA nucleotide
    for nucleotide in sequence:
        if nucleotide not in "ACGT":
            return False

    return True

def update_kmer_count(kmer_data, kmer, next_char):
    """
    Updates the count of a k-mer and tracks the frequency
    of the character that follows it.

    Args:
        kmer_data (dict): Dictionary storing k-mer information
        kmer (str): The current k-mer
        next_char (str): The character that follows the k-mer

    Returns:
        dict: Updated k-mer data dictionary
    """

    # If k-mer is seen for the first time, initialize its data
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}

    # Increment total count for this k-mer
    kmer_data[kmer]['count'] += 1

    # Initialize next character count if not already present
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0

    # Increment frequency of the next character
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k):
    """
    Extracts all k-mers from a sequence and records
    how often each k-mer appears along with the frequency
    of the character that follows it.

    Args:
        sequence (str): DNA sequence
        k (int): Length of k-mers

    Returns:
        dict: Dictionary of k-mer counts and next character frequencies
    """

    kmer_data = {}

    # Loop through sequence, stopping where a k-mer still has a next character
    for i in range(len(sequence) - k):
        # Extract k-mer of length k
        kmer = sequence[i:i+k]

        # Get the character immediately following the k-mer
        next_char = sequence[i+k]

        # Update counts and frequencies
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)

    return kmer_data

def write_results_to_file(kmer_data, output_filename):
    """
    Writes k-mer counts and next character frequencies to a file.

    Output format:
    kmer char1:count char2:count ...

    Args:
        kmer_data (dict): Dictionary of k-mer data
        output_filename (str): File to write results to
    """

    # Sort k-mers alphabetically for consistent output
    sorted_kmers = sorted(kmer_data.keys())

    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            next_chars = kmer_data[kmer]['next_chars']

            # Format next character frequencies
            next_char_str = " ".join(
                f"{char}:{freq}"
                for char, freq in sorted(next_chars.items())
            )

            # Write line to file
            f.write(f"{kmer} {next_char_str}\n")

def main():
    """
    Main function that:
    - Reads input file of DNA sequences
    - Validates each sequence
    - Computes k-mer statistics
    - Writes results to output file

    Command-line arguments:
        1. Input file name
        2. k-mer length
        3. Output file name
    """

    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]

    print(f"Reading sequences from {sequence_file}...")

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            # Skip invalid sequences
            if not validate_sequence(sequence, k):
                print("  Warning: Skipping sequence")
                continue

            # Count k-mers for the sequence
            kmer_data = count_kmers_with_context(sequence, k)

            # Write results to file
            write_results_to_file(kmer_data, output_file)

if __name__ == '__main__':
    main()
