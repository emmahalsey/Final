from kmer_analyzer import validate_sequence, count_kmers_with_context

def test_valid_sequence():
    assert validate_sequence("ATCG", 2) == True

def test_invalid_sequence_length():
    assert validate_sequence("A", 2) == False

def test_invalid_sequence_characters():
    assert validate_sequence("AT1G", 2) == False

def test_kmer_count_basic():
    result = count_kmers_with_context("ATCG", 2)
    assert "AT" in result
