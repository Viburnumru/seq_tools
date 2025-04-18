import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from seq_tools import filter_fastq
import logging
from pathlib import Path


def write_fastq(records, path):
    with open(path, "w") as f:
        SeqIO.write(records, f, "fastq")

def read_fastq(path):
    with open(path) as f:
        return list(SeqIO.parse(f, "fastq"))


def test_filter_by_length(tmp_path):
    record1 = SeqRecord(Seq("ATGC"), id="short", letter_annotations={"phred_quality": [40]*4})
    record2 = SeqRecord(Seq("ATGCGTAT"), id="long", letter_annotations={"phred_quality": [40]*8})
    input_path = tmp_path / "length_test.fastq"
    write_fastq([record1, record2], input_path)

    filter_fastq(str(input_path), "out.fastq", length_bounds=(5, 10))

    records = read_fastq("filtered/out.fastq")
    assert len(records) == 1
    assert records[0].id == "long"

def test_filter_by_quality(tmp_path):
    record1 = SeqRecord(Seq("ATGCGA"), id="lowQ", letter_annotations={"phred_quality": [10]*6})
    record2 = SeqRecord(Seq("GGCCTA"), id="highQ", letter_annotations={"phred_quality": [40]*6})
    input_path = tmp_path / "quality_test.fastq"
    write_fastq([record1, record2], input_path)

    filter_fastq(str(input_path), "out.fastq", quality_threshold=30)

    records = read_fastq("filtered/out.fastq")
    assert len(records) == 1
    assert records[0].id == "highQ"

def test_filter_by_gc_lower_bound(tmp_path):
    record1 = SeqRecord(Seq("ATATAT"), id="lowGC", letter_annotations={"phred_quality": [40]*6})
    record2 = SeqRecord(Seq("GCGCGC"), id="highGC", letter_annotations={"phred_quality": [40]*6})
    input_path = tmp_path / "gc_low_test.fastq"
    write_fastq([record1, record2], input_path)

    filter_fastq(str(input_path), "out.fastq", gc_bounds=(50, 100))

    records = read_fastq("filtered/out.fastq")
    assert len(records) == 1
    assert records[0].id == "highGC"

def test_filter_by_gc_upper_bound(tmp_path):
    record1 = SeqRecord(Seq("GCGCGC"), id="tooGC", letter_annotations={"phred_quality": [40]*6})
    record2 = SeqRecord(Seq("ATGCAT"), id="okGC", letter_annotations={"phred_quality": [40]*6})
    input_path = tmp_path / "gc_upper_test.fastq"
    write_fastq([record1, record2], input_path)

    filter_fastq(str(input_path), "out.fastq", gc_bounds=(0, 70))

    records = read_fastq("filtered/out.fastq")
    assert len(records) == 1
    assert records[0].id == "okGC"


def test_filter_all_pass(tmp_path):
    record = SeqRecord(Seq("ATGCGA"), id="good", letter_annotations={"phred_quality": [40]*6})
    input_path = tmp_path / "pass_test.fastq"
    write_fastq([record], input_path)

    filter_fastq(str(input_path), "out.fastq", length_bounds=(1,10), quality_threshold=30, gc_bounds=(0,100))

    records = read_fastq("filtered/out.fastq")
    assert len(records) == 1
    assert records[0].id == "good"
    

def test_empty_input_file(tmp_path):
    input_path = tmp_path / "empty.fastq"
    input_path.write_text("") 
    output_file = "empty_out.fastq"

    filter_fastq(str(input_path), output_file)

    output_path = Path("filtered") / output_file
    assert not output_path.exists() 


def test_no_output_file_if_no_records_pass(tmp_path):
    record = SeqRecord(Seq("ATATAT"), id="low_quality", letter_annotations={"phred_quality": [5] * 6})
    input_path = tmp_path / "input.fastq"
    write_fastq([record], input_path)

    output_file = "no_output.fastq"
    filter_fastq(str(input_path), output_file, quality_threshold=40, gc_bounds=(80, 100)) 

    output_path = Path("filtered") / output_file
    assert not output_path.exists(), "Файл не должен быть создан, если ни одна запись не прошла фильтрацию"


def test_invalid_file_format(tmp_path):
    input_file = tmp_path / "not_a_fastq.txt"
    input_file.write_text("This is not a fastq file")

    with pytest.raises(ValueError, match="Input file must be in FASTQ format"):
        filter_fastq(str(input_file), "out.fastq")