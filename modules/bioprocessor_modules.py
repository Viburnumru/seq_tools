def genes_from_gbk(path_gbk: str) -> list[dict[str:str]]:
    """extracts gene names and their translations from a GenBank (.gbk) file.
    This function reads a GenBank file and
    collects genes and their corresponding
    protein translations from CDS regions.
    It stores each gene and its translation
    in a dictionary, which is then added to a list.
    Args:
        path_gbk (str): path to input gbk file

    Returns:
        List[Dict[str, str]]: A list of dictionaries
        where each dictionary contains
        a gene name ('gene') and its translation ('translation').
    """
    genes_info = []
    current_gene = None

    with open(path_gbk, "r") as file:
        for line in file:
            if line.strip().startswith("CDS"):
                if current_gene and current_gene["gene"] is not None:
                    genes_info.append(current_gene)
                current_gene = {"gene": None, "translation": ""}
            elif "/gene" in line and current_gene:
                gene_name = line.split("=")[1].replace('"', "").strip()
                current_gene["gene"] = gene_name
            elif (
                "/translation=" in line
                and current_gene
                and current_gene["gene"] is not None
            ):
                translation = line.split("=")[1].replace('"', "").strip()
                current_gene["translation"] = translation
            elif (
                not line.strip().startswith("/")
                and current_gene
                and current_gene["translation"] is not None
            ):
                current_gene["translation"] += line.strip().replace('"', "")
        if current_gene:
            genes_info.append(current_gene)
    return genes_info


def find_genes_of_interest(
    genes_info: list[dict[str:str]],
    genes_of_interest: list,
    n_before: int,
    n_after: int,
) -> list[dict[str:str]]:
    """_summary_

    Args:
        genes_info (list[dict[str:str]]): A list of dictionaries
        where each dictionary contains
        a gene name ('gene') and its translation ('translation').
        genes_of_interest (list): list with genes names to search
        n_before (int): number of genes to search before gene of interest
        n_after (int): number of genes to search after gene of interest

    Returns:
        list[dict[str:str]]: A list of dictionaries
        where each dictionary contains
        a neighbor gene name ('gene') and its translation ('translation').
    """
    found_genes = []
    for idx, gene in enumerate(genes_info):
        gene_name = gene["gene"]
        if gene_name is not None and any(
            gene_of_interest in gene_name
            for gene_of_interest in genes_of_interest
        ):
            if idx == 0:  # first gene
                start_idx = 0
                end_idx = min(len(genes_info), idx + n_after + 1)
                for neighbor_idx in range(idx + 1, end_idx):
                    if genes_info[neighbor_idx] not in found_genes:
                        found_genes.append(genes_info[neighbor_idx])
            elif idx == len(genes_info) - 1:  # last gene
                start_idx = max(0, idx - n_before)
                end_idx = len(genes_info)
                for neighbor_idx in range(start_idx, idx):
                    if genes_info[neighbor_idx] not in found_genes:
                        found_genes.append(genes_info[neighbor_idx])
            else:  # middle gene
                start_idx = max(0, idx - n_before)
                end_idx = min(len(genes_info), idx + n_after + 1)
                for neighbor_idx in range(start_idx, end_idx):
                    if (
                        neighbor_idx != idx
                        and genes_info[neighbor_idx] not in found_genes
                    ):
                        found_genes.append(genes_info[neighbor_idx])

    return found_genes


def save_to_fasta(genes: list[dict[str:str]], output_fasta: str) -> None:
    """save found genes from list of dictionaries to fasta format.
    Results are saved in output file.

    Args:
        genes (list[dict[str:str]]): A list of dictionaries
        where each dictionary contains
        a gene name ('gene') and its translation ('translation').
        output_fasta (str): path to output file
    """
    with open(output_fasta, "w") as file:
        for gene in genes:
            if gene["translation"]:
                file.write(f">{gene['gene']}\n")
                file.write(f"{gene['translation']}\n")
