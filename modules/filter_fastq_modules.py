import numbers
import os


def make_bounds(bounds: int, tuple) -> tuple:
    """Converts a given bounds into a tuple format.

    Args:
        bounds (int, tuple): given upper bound or
        lower and upper bounds

    Returns:
        tuple: bounds (lower, upper)
        in case lower bound was not given returns (0, upper)
    """
    if isinstance(bounds, numbers.Real):
        return (0, bounds)
    return bounds


def is_bounded(x: int, bounds: tuple) -> bool:
    """Checks if number lays between given bounds

    Args:
        x (int): number of interest
        bounds (tuple): bounds (lower, upper)

    Returns:
        bool: True if number within bounds
    """
    return bounds[0] <= x <= bounds[1]


def gc_count(na: str) -> float:
    """
    Counts GC in seq.
    Args:
        na (string): NA sequence
    Returns:
        float: (%) GC content
    """
    na = na.upper()
    return 100 * (na.count("G") + na.count("C")) / len(na)


def is_gc_valid(
    sequence: str, gc_bounds: tuple[float, float] = (0, 100)
) -> bool:
    """Checks if sequence valid according GC content

    Args:
        sequence (str): sequence of interest
        gc_bounds (tuple[float, float], optional):
        GC bounds. Defaults to (0, 100).

    Returns:
        bool: True if GC content of sequence is within GC bounds
    """
    gc_bounds = make_bounds(gc_bounds)
    gc = gc_count(sequence)

    return is_bounded(gc, gc_bounds)


def is_length_valid(
    sequence: str, length_bounds: tuple[float, float] = (0, 2**32)
) -> bool:
    """Checks if sequence valid according sequence length

    Args:
        sequence (str): sequence of interest
        length_bounds (tuple[float, float], optional):
        length bounds. Defaults to (0, 2**32).

    Returns:
        bool: True if sequence length is within length bounds
    """

    length_bounds = make_bounds(length_bounds)
    return is_bounded(len(sequence), length_bounds)


def is_quality_valid(quality: str, quality_threshold: float = 0) -> bool:
    """Checking quality of the sequence

    Args:
        quality (str): sequence of interest
        quality_threshold (float, optional): quality bounds. Defaults to 0.

    Returns:
        bool: True if sequence quality is within quality bounds
    """

    average_score = sum(ord(char) - 33 for char in quality) / len(quality)
    return average_score >= quality_threshold


def read_fastq_file(path: str) -> dict[str, tuple[str, str]]:
    """read FASTQ file and converts data into dictionary

    Args:
        path (str): path to file

    Returns:
        dict[str, tuple[str, str]]: dictionary[name, (sequence, quality)]
    """
    with open(path, "r") as file:
        seqs = {}
        while True:
            line = file.readline().strip()
            if not line:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()
            name = line[0:].split(" ")[0]
            seqs[name] = (sequence, quality)
    return seqs


def write_output_fastq(
    name: str, sequence: str, quality: str, output_file: str
) -> None:
    """Append a single sequence to the output FASTQ file
    in the 'filtered' directory.
    Args:
        name (str): the sequence identifier (e.g., @SEQ_ID in FASTQ format).
        sequence (str): sequence
        quality (str): quality in FASTQ format
        output_path (str): path to output file
    """
    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    with open(output_path, "a") as file:
        file.write(f"{name}\n")
        file.write(f"{sequence}\n")
        file.write("+\n")
        file.write(f"{quality}\n")
