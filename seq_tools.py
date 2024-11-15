from modules.filter_fastq_modules import (
    is_gc_valid,
    is_length_valid,
    is_quality_valid,
    read_fastq_file,
    write_output_fastq,
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
    path: str,
    output_file: str = "output_filter_fastq.fastq",
    gc_bounds: tuple[float, float] = (0, 100),
    length_bounds: tuple[float, float] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """filter FASTQ sequences according to length, GC and quality thresholds.
    Result is saved in output file in a direcroty 'filtered'.

    Args:
        path (str): path to FASTQ sequences
        output_file (str, optional): output file name with extension.
        Defaults to 'output_filter_fastq.fastq'.
        gc_bounds (tuple[float, float], optional):
        GC bounds. Defaults to (0, 100).
        length_bounds (tuple[float, float], optional):
        length bounds. Defaults to (0, 2**32).
        quality_threshold (float, optional): quality bound. Defaults to 0.
    """
    seqs = read_fastq_file(path)

    for name, (sequence, quality) in seqs.items():
        if (
            is_gc_valid(sequence, gc_bounds)
            and is_length_valid(sequence, length_bounds)
            and is_quality_valid(quality, quality_threshold)
        ):
            write_output_fastq(name, sequence, quality, output_file)


def run_dna_rna_tools(*args: str) -> list:
    """finding different chains of DNA and/or RNA.
    In case of transmission incorrect type of operation or RNA/DNA sequence,
    the corresponding error will be printed on the screen.

    Args:
        positional arguments, consisting of nucleic acid sequences.
        Last argument is operation:
        transcribe, reverse, complement or reverse complement

    Returns:
        list: list of strings sequences or None's (in case of errors)
        str: if there is only 1 result

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
            ans.append(None)
            stderr[seq] = (0, "operation type is not supported")
            continue
        if not is_na(seq):
            ans.append(None)
            stderr[seq] = (0, "is not NA")
            continue
        if not na_type(seq):
            ans.append(None)
            stderr[seq] = (na_type(seq), "RNA and DNA mix")
        else:
            if task != "transcribe":
                ans.append(operations.get(task)(seq))
            elif task == "transcribe" and na_type(seq) == "DNA":
                ans.append(operations.get(task)(seq))
            else:
                ans.append(None)
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
