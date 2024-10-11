## Seq_tools

Seq Tools is a toolkit for working with DNA and RNA sequences, as well as for filtering FASTQ data. The project includes functions for finding the complementary, reverse, and reverse complementary strands of DNA or RNA, transcribing DNA, and filtering sequences based on various criteria.

 ```python
    -/
     |- README.md
     |- seq_tools.py
     |- bio_files_processor.py
     |- modules/
           |- dna_rna_tools_modules.py
           |- filter_fastq_modules.py
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
To run `Seq_tools` or `bio_files_processor.py` you need to have Python 3.x installed. To check this: 
 ```bash
python --version
  ```
clone or download manually repository: 
 ```
git clone  https://github.com/Viburnumru/seq_tools.git
 ```

## Main functions
### Seq_tools

Main script `seq_tools.py` consist of `filter_fastq` and `run_dna_rna_tools` functions.

`Run_dna_rna_tools` performs the specified operations on a DNA or/and RNA sequence.  
Operations includes "transcribe", "reverse", "complement" and "reverse_complement".
In case of attempting to pass an incorrect operation type or RNA/DNA sequence, an appropriate error will be printed on the screen.

Returns:
list: a list of results  
str: a single result  

`Filter_fastq` filters reads file in FASTQ format and filters sequences based on GC content, read length, and read quality.

Filtered sequences are stored in a file in the folder 'filtered'. Filtered sequences are the same structure as  input FASTQ, but only sequences passed all criteria are included.

### Bio_files_processor

Main script `Bio_files_processor.py` operates with Fasta and Gbk formats. It allows conversion of multiline FASTA sequences to one row, selects genes from gbk and stores them in fasta formst, and parses blast results to find out best aligned proteins 


## Modules

Additional functions are stored separately in the folder 'modules':
1. **dna_rna_tools_modules.py**
   consist of functions required for sequence check:
    - *is_na* checks if sequence contains nucleotides;
    - *na_type* checks if sequence contains only one type (DNA or RNA) nucleotides;
      
   and main functions:
   
    - funcion *transcribe* (works only with DNA sequences);
    - function *complement* returns complementary strand; 
    - function *reverse returns* reverse strand;
    - function *reverse_complement* returns reverse complement strand;
 These functions are case sensitive.

2. **filter_fastq_modules.py**
   consist of functions required for boundaries check:
   - *make_bounds* make bounds in required tuple format
   - *is_bounded* checks if criteria passed or not
  
     
   and mail functions:   
   - *gc_count*, which returns GC content in %;
   - *filter_gc*, which checks if sequence passes gc threshold;
   - *filter_length*, which checks if sequence passes length threshold;
   - *filter_quality*, which checks if sequence passes quality threshold

These functions requires FASTQ input format (see *Examples of usage*).

3. **bioprocessor_modules.py**
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
filter_fastq(path, 'example_output.fastq', gc_bounds=(0, 100), length_bounds=(0,2**32), quality_threshold=0)
```

**Seqtools: dna_rna_tools**  
Input:  
```
run_dna_rna_tools('AtgC', 'AUt', 'Anna', 'gaC', 'reverse')
```
Output:  
```
['CgtA', None, None, 'Cag']
```
on the screen following errors were printed:  
```
 AUt : RNA and DNA mix
 Anna : is not NA
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
 
