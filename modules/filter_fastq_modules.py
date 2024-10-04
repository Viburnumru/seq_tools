def gc_count(na: str) -> float:
    """
    Функция считает GC-состав.
    Параметры: строка 'na'.
    Возвращает значение с float (в процентах).
    """
    na = na.upper()
    return 100 * (na.count("G") + na.count("C")) / len(na)


def filter_gc(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: tuple[float, float] = (0, 100)
) -> dict[str, tuple[str, str]]:
    """
    Функция фильтрует последовательность нуклеотидов по GC-составу.
    Параметры:
    seqs: dict[str, tuple[str, str]] -
    словарь, ключами которого являются названия ридов,
    а значениями - кортежи, состоящие из сиквенса и его качества.
    gc_bounds (tuple[float, float]): Интервал GC-состава для фильтрации
    (по умолчанию (0, 100)).
    Может приниматься только верхняя граница в виде float.
    Тогда по умолчанюю нижняя граница считается равной 0.
    Возвращает словарь вида dict[str, tuple[str, str]],
    ключами которого являются названия сиквенсов,
    отфильтрованных по GC-составу ридов,
    а значениями - кортежи, состоящие из сиквенса и его качества.

    """
    filtered_seqs = {}
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    for name, (sequence, quality) in seqs.items():
        gc = gc_count(sequence)
        if gc_bounds[0] <= gc <= gc_bounds[1]:
            filtered_seqs[name] = (sequence, quality)
    return filtered_seqs


def filter_length(
    seqs: dict[str, tuple[str, str]],
    length_bounds: tuple[float, float] = (0, 2**32)
) -> dict[str, tuple[str, str]]:
    """
    Функция фильтрует последовательность нуклеотидов по длине рида.
    Параметры:
    seqs: dict[str, tuple[str, str]] - словарь, где ключи - это названия ридов,
    а значения - кортежи, состоящие из прочтения и его качества.
    length_bounds (tuple[int, int]): Интервал длины для фильтрации
    (по умолчанию (0, 2**32)).
    Может приниматься только верхняя граница в виде float.
    Тогда по умолчанюю нижняя граница считается равной 0.
    Возвращает словарь вида dict[str, tuple[str, str]],
    ключами которого являются названия сиквенсов,
    отфильтрованных по длине рида,
    а значениями - кортежи, состоящие из сиквенса и его качества.

    """
    filtered_seqs = {}
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)
    for name, (sequence, quality) in seqs.items():
        if length_bounds[0] <= len(sequence) <= length_bounds[1]:
            filtered_seqs[name] = (sequence, quality)
    return filtered_seqs


def filter_quality(
    seqs: dict[str, tuple[str, str]], quality_threshold: float = 0
) -> dict[str, tuple[str, str]]:
    """
    Функция фильтрует последовательность нуклеотидов
    по заданному уровню качества рида.
    Параметры:
    seqs: dict[str, tuple[str, str]] - словарь,
    ключами которого являются названия ридов,
    а значениями - кортежи, состоящие из сиквенса и его качества.
    quality_threshold (float): Пороговое значение среднего качества рида,
    по умолчанию 0.
    Возвращает словарь вида dict[str, tuple[str, str]],
    ключами которого являются названия сиквенсов, отфильтрованных по качеству,
    а значениями - кортежи, состоящие из сиквенса и его качества.

    """
    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():
        quality_scores = []
        for char in quality:
            quality_scores.append(ord(char) - 33)
        average_score = sum(quality_scores) / len(quality_scores)
        if average_score >= quality_threshold:
            filtered_seqs[name] = (sequence, quality)
    return filtered_seqs
