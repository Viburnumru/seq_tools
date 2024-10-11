from modules.bioprocessor_modules import (
    genes_from_gbk,
    find_genes_of_interest,
    save_to_fasta,
)


def convert_multiline_fasta_to_oneline(
    input_file: str, output_file: str = "output_fasta_conversion.fastq"
) -> None:
    """converts multiline FASTA sequence to one line format

    Args:
        input_file (str): path to input file
        output_file (str, optional): path to output file.
        Defaults to 'output_fasta_conversion.fastq'.
    """
    with open(input_file, "r") as input, open(output_file, "w") as output:
        name = None
        sequence = []
        for line in input:
            line = line.strip()
            if line.startswith(">"):

                if name:
                    output.write(name + "\n")
                    output.write("".join(sequence) + "\n")

                name = line
                sequence = []
            else:
                sequence.append(line)

        if name:
            output.write(name + "\n")
            output.write("".join(sequence) + "\n")


def parse_blast_output(
    input_file: str, output_file: str = "output_blast.txt"
) -> None:
    """parses the names of the best proteins
    from the txt file of the blast results.
    Proteins are seved in output file.

    Args:
        input_file (str): path to input txt file eith blast results
        output_file (str, optional): path to output file.
        Defaults to 'output_blast.txt'.
    """
    with open(input_file, "r") as file:
        unique_descriptions = set()
        current = False
        for line in file:
            line = line.strip()

            if line.startswith("Description"):
                current = True
                continue

            if current and line.startswith(">"):
                description = line[1:].split("[")[0].strip()
                if "MULTISPECIES: " in description:
                    description = description.replace("MULTISPECIES: ", "")
                if description not in unique_descriptions:
                    unique_descriptions.add(description)
                    current = False

        sorted_descriptions = sorted(unique_descriptions, key=str.lower)

    with open(output_file, "w") as file:
        for description in sorted_descriptions:
            file.write(f"{description}\n")


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: list,
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "neighbor_genes_from_gbk_search.fasta",
) -> None:
    """Select neigbors of genes of interest from a GenBank
    (.gbk) file and save them into a file.
    Args:
        input_gbk (str): path to input (.gbk) file
        genes (list): list of genes names to search
        n_before (int): number of genes to search before gene of interest
        n_after (int): number of genes to search after gene of interest
        output_fasta (str): path to output file.
        Defaults to "neighbor_genes_from_gbk_search.fasta".

    """
    genes_info = genes_from_gbk(input_gbk)
    found_genes = find_genes_of_interest(genes_info, genes, n_before, n_after)
    save_to_fasta(found_genes, output_fasta)
