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


# 🔽 ADD NEW TESTS BELOW THIS LINE

def test_kmer_counts_exact():
    result = count_kmers_with_context("ATAT", 2)
    assert result["AT"]["count"] == 1

def test_empty_sequence():
    assert count_kmers_with_context("", 2) == {}

def test_k_equals_length():
    result = count_kmers_with_context("ATCG", 4)
    assert result == {}

def test_kmer_full_output():
    result = count_kmers_with_context("ATCG", 2)
    
    # Expected:
    # AT -> C
    # TC -> G
    
    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["C"] == 1
    
    assert result["TC"]["count"] == 1
    assert result["TC"]["next_chars"]["G"] == 1
