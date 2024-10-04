def is_na(seq_na: str) -> bool:
    """
    Функция проверяет,
    является ли строка нуклеиновой кислотой.

    Args:
        seq_na (str): последовательность

    Returns:
        bool: True, если последовательность
        состоит только из нуклеотидов, иначе False.
    """
    na_list = ["A", "T", "G", "C", "U", "a", "t", "g", "c", "u", "a"]
    return all(ch in na_list for ch in seq_na)


def na_type(na: str) -> str or bool:
    """
    Функция определяет тип нуклеиновой кислоты.

    Args:
        na (str): последовательность

    Returns:
        str or bool: возвращает False, если тип нельзя определить;
        возвращает строку 'RNA', если последовательность содержит 'U',
        'DNA' в остальных случаях (по умолчанию).
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
    Функция возвращает транскрипт цепи ДНК

    Args:
        na (str): последовательность ДНК

    Returns:
        str: последовательность РНК
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
    Функция находит обратную цепь ДНК или РНК.

    Args:
        na (str): последовательность.

    Returns:
        str: обратная последовательность.
    """
    return na[::-1]


def complement(na: str) -> str:
    """
    Функция находит комплементарную цепь ДНК или РНК.

    Args:
        na (str): последовательность.

    Returns:
        str: комплементарная последовательность.
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
    Функция находит обратную комплементарную цепь ДНК или РНК.

    Args:
        na (str): последовательность.

    Returns:
        str: обратная комплементарная последовательность.
    """
    return reverse(complement(na))
