from typing import Union


def is_na(seq_na: str) -> bool:
    """
    Checking if string is nucleic acid (NA)

    Args:
        seq_na (str): NA sequence

    Returns:
        bool: True, if sequence consist from nucleotides only
    """
    na_list = ["A", "T", "G", "C", "U", "a", "t", "g", "c", "u", "a"]
    return all(ch in na_list for ch in seq_na)


def na_type(na: str) -> Union[str, bool]:
    """
    Finding of NA type

    Args:
        na (str): sequence

    Returns:
        str or bool: returns False, if NA type could not be defined;
        string 'RNA', if sequence contain 'U',
        otherwise 'DNA' (default).
    """
    for base in na:
        has_t = "T" in na.upper()
        has_u = "U" in na.upper()
        if has_t and has_u:
            return False
        if has_t:
            return "DNA"
        return "RNA"


def transcribe(na: str) -> str:
    """
    Finding DNA trancscripte

    Args:
        na (str): DNA sequence

    Returns:
        str: RNA sequence
    """
    dna_to_rna_dict = {
        "A": "A",
        "T": "U",
        "G": "G",
        "C": "C",
        "a": "a",
        "t": "u",
        "g": "g",
        "c": "c",
    }
    return "".join(dna_to_rna_dict[base] for base in na)


def reverse(na: str) -> str:
    """
    Finding reverse chain of DNA or RNA

    Args:
        na (str): sequence

    Returns:
        str: reverse sequence
    """
    return na[::-1]


def complement(na: str) -> str:
    """
     Finding complementary chain of DNA or RNA

    Args:
        na (str): sequence

    Returns:
        str: complementary sequence
    """
    na_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "U": "A",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "u": "a",
        "c": "g",
    }
    return "".join(na_dict[base] for base in na)


def reverse_complement(na: str) -> str:
    """
    Finding reverse complement of DNA or RNA

    Args:
        na (str): sequence

    Returns:
        str: reverse complement
    """
    return reverse(complement(na))
