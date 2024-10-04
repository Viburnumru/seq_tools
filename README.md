## Seq_tools

Seq Tools is a toolkit for working with DNA and RNA sequences, as well as for filtering FASTQ data. The project includes functions for finding the complementary, reverse, and reverse complementary strands of DNA or RNA, transcribing DNA, and filtering sequences based on various criteria.

 ```python
    -/
     |- README.md
     |- seq_tools.py 
     |- modules/
           |- dna_rna_tools_modules.py
           |- filter_fastq_modules.py
 ```

Author: Anna Kalinina

## Contents

- [Installation](#Installation)
- [Usage](##Usage)
  - [Main functions](###Main-functions)
- [Modules](##modules)
- [Examples](##examples)

## Installation
To run `Seq_tools` you need to have Python 3.x installed. To check this: 
 ```bash
python --version
  ```
clone or download manually repository: 
 ```
git clone  https://github.com/Viburnumru/seq_tools.git
 ```
The script is not supposed to be used from commandline. Use any Python IDEs to run `seq_tools`.
You may work directly in the file `seq_tools.py` or import needed modules into your script.

## Usage

### Main functions

Main script consist of `filter_fastq` and `run_dna_rna_tools` functions.

`Run_dna_rna_tools` performs the specified operations on a DNA or\and RNA sequence.
args: positional arguments, strings representing DNA/RNA sequences, the last argument must be a string indicating the type of operation (for example, "transcribe", "reverse", etc.). In case of attempting to pass an incorrect operation type or RNA/DNA sequence, an appropriate error will be printed on the screen.

Returns:
list: a list of results
str: a single result

`Filter_fastq` filters reads in FASTQ format based on GC content, read length, and read quality.
Parameters:
seqs: a dictionary where the keys are read names, and the values are tuples consisting of the sequence and its quality.
gc_bounds: the GC content interval for filtering (default is (0, 100)).
length_bounds: the length interval for filtering (default is (0, 2**32)).
quality_threshold (float): the threshold value for the average quality of the read, default is 0.
Returns:
A dictionary with the same structure as seqs.

## Modules

Additional functions are stored separately in the folder 'modules':
1. dna_rna_tools_modules.py
   consist of functions required for sequence check:
    - is_na checks if sequence contains nucleotides;
    - na_type checks if sequence contains only one type (DNA or RNA) nucleotides;
      
   and main functions:
   
    - funcion transcribe (works only with DNA sequences);
    - function complement;
    - function reverse;
    - function reverse complement;

2. filter_fastq_modules.py
   consist of
   - gc_count, which returns GC content in %
   - filter_gc
   - filter_length
   - filter_quality


All functions from modules could be imported and used separately.


## Examples
```

EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079805': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
    '@SRX079806': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT', 'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
    '@SRX079807': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT', 'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
    '@SRX079808': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', 'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
    '@SRX079809': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
    '@SRX079810': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
    '@SRX079811': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA', '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
    '@SRX079812': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG', '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
    }  
```

Input: 

```
filter_fastq(FASTQ, gc_bounds=(40,90), length_bounds=(40,80), quality_threshold=33)
```
Output: 
```
{'@SRX079806': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT',
  'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
 '@SRX079810': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC',
  'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE')}
```


Input:
```
run_dna_rna_tools('AtgC', 'AUt', 'Anna', 'gaC', 'reverse')
```
Output:
```
['CgtA', 'Cag']
```
on the screen following errors were printed:

AUt : RNA and DNA mix  
Anna : is not NA




