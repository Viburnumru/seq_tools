## Seq_tools

Seq Tools is a toolkit for working with DNA and RNA sequences, as well as for filtering FASTQ data. The project includes functions for finding the complementary, reverse, and reverse complementary strands of DNA or RNA, transcribing DNA, and filtering sequences based on various criteria.

 ```python
    -/
     |- README.md
     |- seq_tools.py
     |- bio_files_processor.py
     |- modules/
           |- bioprocessor_modules.py
     |- examples/
 ```

Author: Anna Kalinina

## Contents

- [Installation](#Installation)
- [Main functions](##Main-functions)
- [Seq tools](###Seq_tools)
- [Bio_files_processor](###Bio_files_processor)
- [Modules](##modules)
- [Examples of usage](##examples-of-usage)

## Installation
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

`Run_dna_rna_tools` performs the specified operations on a DNA or/and RNA sequence.
Operations includes "transcribe", "reverse", "complement" and "reverse_complement". In case of attempting to pass an incorrect operation type or RNA/DNA sequence, an appropriate error will be printed on the screen.

Returns: list: a list of results
str: a single result

`Filter_fastq` filters reads file in FASTQ format and filters sequences based on GC content, read length, and read quality.

Filtered sequences are stored in a file in the folder 'filtered'. Filtered sequences are the same structure as  input FASTQ, but only sequences passed all criteria are included. If the filtered directory does not exist, it will be created automatically.

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

**Seqtools: dna_rna_tools**  
Input:  
```
DNASequence("atCG").complement()
```
Output:  taGC

The acces to previous functionallity is still available. In case you need to manipulate several sequences simultaneously: 
Input:
```
run_dna_rna_tools('AtgC', 'AUt', 'Anna', 'gaC', 'reverse_complement')

```
Output:
```
['GcaT', None, None, 'Gtc']
```
Following errors will be printed on screen:
```
AUt : is not a valid DNA or RNA sequence
Anna : is not a valid DNA or RNA sequence
```

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
 
