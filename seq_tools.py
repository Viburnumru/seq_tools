from modules.filter_fastq_modules import (
    filter_gc,
    filter_length,
    filter_quality,
)
from modules.dna_rna_tools_modules import (
    is_na,
    na_type,
    transcribe,
    reverse,
    complement,
    reverse_complement,
)


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: tuple[float, float] = (0, 100),
    length_bounds: tuple[float, float] = (0, 2**32),
    quality_threshold: float = 0,
) -> dict[str, tuple[str, str]]:
    """
    Фильтрует риды в формате FASTQ по GC-составу, длине рида и качеству рида.
    Параметры:
    seqs (dict[str, tuple[str, str]]): словарь,
    ключами которого являются названия ридов,
    а значениями - кортежи, состоящие из сиквенса и его качества.
    gc_bounds (tuple[float, float]): Интервал GC-состава для фильтрации
    (по умолчанию (0, 100)).
    length_bounds (tuple[int, int]): Интервал длины для фильтрации
    (по умолчанию (0, 2**32)).
    quality_threshold (float): Пороговое значение среднего качества рида,
    по умолчанию 0.
    Возвращает словарь вида dict[str, tuple[str, str]],
    ключами которого являются названия отфильтрованных сиквенсов,
    а значениями - кортежи, состоящие из сиквенса и его качества.

    """

    filtered_seqs = filter_gc(seqs, gc_bounds)
    filtered_seqs = filter_length(filtered_seqs, length_bounds)
    filtered_seqs = filter_quality(filtered_seqs, quality_threshold)
    return filtered_seqs


def run_dna_rna_tools(*args: str) -> list:
    """
    Функция выполняет заданные операции с последовательностью ДНК или РНК.
    Параметры (args):
        - позиционные аргументы: строки,
        представляющие собой последовательности ДНК/РНК,
        - последний аргумент должен быть строкой,
        указывающей на тип операции
        (например, "transcribe", "reverse", и т.д.).
        В случае попытки передать неправильный
        тип операции или последовательности РНК/ДНК,
        соответствующая ошибка будет напечатана на экране.

    Returns:
        list: список результатов
        str: только один результат
    """
    task = args[-1]
    operations = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    stderr = {}
    ans = []
    for seq in args[:-1]:
        if task not in operations:
            stderr[seq] = (0, "operation type is not supported")
        else:
            if not is_na(seq):
                stderr[seq] = (is_na(seq), "is not NA")
            else:
                if not na_type(seq):
                    stderr[seq] = (na_type(seq), "RNA and DNA mix")
                else:
                    if task != "transcribe":
                        ans.append(operations.get(task)(seq))
                    elif task == "transcribe" and na_type(seq) == "DNA":
                        ans.append(operations.get(task)(seq))
                    else:
                        stderr[seq] = (
                            0,
                            "operation type is not compatible with NA type",
                        )
    if stderr:
        for seq, (result, message) in stderr.items():
            print(f" {seq} : {message}")
    if len(ans) == 1:
        return ans[0]
    return ans
