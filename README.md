## Seq_tools

Seq Tools is a toolkit for working with DNA and RNA sequences, as well as for filtering FASTQ data. The project includes functions for finding the complementary, reverse, and reverse complementary strands of DNA or RNA, transcribing DNA, and filtering sequences based on various criteria.

 ```python
    -/
     |- README.md
     |- seq_tools.py
     |- bio_files_processor.py
     |- fastq_filter.log
     |- requirements.txt
     |- modules/
           |- bioprocessor_modules.py
     |- examples/
     |- tests/
            | - test_fastq_filter.py

 ```

Author: Anna Kalinina

## Contents

- [Installation](#Installation)
- [Main functions](##Main-functions)
- [Seq tools](###Seq_tools)
- [Bio_files_processor](###Bio_files_processor)
- [Modules](##modules)
- [Examples of usage](##examples-of-usage)
- [Examples of usage from command line](###examples-of-usage-from-command-line)
- [Tests](##tests)
## Installation and requirements
To run `Seq_tools` or `bio_files_processor.py` you need to have Python 3.x, biopython (Bio) module installed. 

clone or download manually repository: 
 ```
git clone  https://github.com/Viburnumru/seq_tools.git
 ```

## Main functions
### Seq_tools

Main script `seq_tools.py`.

`seq_tools.py` performs the specified operations on a DNA or/and RNA sequence.  
Operations includes "transcribe" for DNA, "reverse", "complement" and "reverse_complement" for both DNA and RNA.

This Python library provides classes for handling and manipulating biological sequences, including DNA, RNA, and protein sequences. It is built using object-oriented programming principles and leverages the abc module for abstract base classes.

Abstract Base Class: BiologicalSequence provides a foundation for all sequence types.

DNA and RNA Sequences: Classes for handling nucleic acid sequences, including complement, reverse, and reverse-complement operations.

Protein Sequences: Class for handling amino acid sequences, including validation and alphabet checking.

`Filter_fastq` filters reads file in FASTQ format and filters sequences based on GC content, read length, and read quality.

Filtered sequences are stored in a file in the folder 'filtered'. Filtered sequences are the same structure as  input FASTQ, but only sequences passed all criteria are included. If the filtered directory does not exist, it will be created automatically.

Filter_fastq is available from command line, see examples of usage below.

**Parameters**:
input_file (str): Path to the input FASTQ file.

output_file (str, optional): Name of the output file. Defaults to output_filter_fastq.fastq.

gc_bounds (tuple[float, float], optional): Minimum and maximum GC content (as percentages). Defaults to (0, 100).

length_bounds (tuple[int, int], optional): Minimum and maximum sequence length. Defaults to (0, 2**32).

quality_threshold (float, optional): Minimum average quality score (Phred score). Defaults to 0.

### Bio_files_processor

Main script `Bio_files_processor.py` operates with Fasta and Gbk formats. It allows conversion of multiline FASTA sequences to one row, selects genes from gbk and stores them in fasta formst, and parses blast results to find out best aligned proteins 


## Modules

Additional functions are stored separately in the folder 'modules':
1.  **bioprocessor_modules.py**
   - *genes_from_gbk* extracts gene names and their translations from a GenBank (.gbk) file.  
   - *find_genes_of_interest* searchs neigbor genes of a target gene  
   - *save_to_fasta* saves genes to fasta file  
   
All functions from modules could be imported and used separately.  

## Examples of usage

You may find example input and output files in the **Examples** folder.  

**Seq_tools: filter_fastq**  
Filtering with default parameters.    
Example file - *example_fastq.fastq*.   
Input:  
```
filter_fastq(path, 'filtered_example.fastq', gc_bounds=(0, 100), length_bounds=(0,2**32), quality_threshold=0)
```
### Example of usage filter fastq_seq from command line

```
python3 seq_tools.py examples/example_fastq.fastq -o filtered_results.fastq --gc 30 70 --length 50 150 --quality 20
```
2025-04-13 01:23:47,206 - INFO - Starting processing of file: examples/example_fastq.fastq


2025-04-13 01:23:47,226 - INFO - Successfully saved 34 sequences to filtered/filtered_results.fastq

or
```
python3 seq_tools.py examples/example_fastq.fastq -o filtered_results1.fastq --gc 90 100 --length 50 150 --quality 100
```
2025-04-15 00:31:47,312 - INFO - Starting processing of file: examples/example_fastq.fastq


2025-04-15 00:31:47,339 - WARNING - No sequences passed the filtering criteria

Logging info is saved in fastq_filter.log. 


**Seqtools:**  
Input:  
```
DNASequence("atCG").complement()
```
Output:  taGC
 


**Bio_files_processor**  

convert_multiline_fasta_to_oneline  
example file - *example_multiline_fasta.fasta  

Input:  
```
convert_multiline_fasta_to_oneline(path_convertion, 'output_file_conversion.fastq')
``` 
example file - *example_gbk.gbk*  

Input:
```
parse_blast_output(path_blast, 'blast_output.txt')
```

select_genes_from_gbk_to_fasta  
example file - *example_gbk.gbk*  

Input:
```
select_genes_from_gbk_to_fasta(path_gbk, genes=['dtpD', 'pxpB'], n_before=2, n_after=2, output_fasta="output_genes_from_gbk.fasta")

```
 
## Tests
Run pytest from ./ folder:
```
pytest
```

seq_tools % pytest                                         
=============================================================================== test session starts ===============================================================================
platform darwin -- Python 3.13.2, pytest-8.3.5, pluggy-1.5.0


rootdir: /#########/seq_tools


collected 8 items                                                                                                                                                                 

tests/test_fastq_filter.py ........                                                                                                                                         [100%]

================================================================================ 8 passed in 0.21s ================================================================================
