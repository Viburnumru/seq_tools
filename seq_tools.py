from abc import ABC, abstractmethod
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """Abstract base class representing a biological sequence (e.g., DNA, RNA, protein).
    This class provides basic functionality for handling biological sequences,
    including validation, length calculation, and indexing. Subclasses must
    implement the `_check_alphabet` method to ensure the sequence contains
    valid characters for the specific type of biological sequence.

    Args:
        sequence (str): The biological sequence as a string.
    Raises:
        ValueError: If the sequence contains invalid characters (as determined
                    by the `_check_alphabet` method).
    """
    def __init__(self, sequence: str):
        """initialization

        Args:
            sequence (str): biosequence as a string

        Raises:
            ValueError: in case of wrong characteres 
        """
        self.sequence = sequence
        if not self._check_alphabet():
            raise ValueError(
                f"Invalid characters in sequence: {self.sequence}"
            )

    def __len__(self) -> int:
        """
        Returns:
            int: returns len of sequence
        """
        return len(self.sequence)

    def __getitem__(self, index):
        """Get a character or slice of the biological sequence by index.

        Args:
            index (int or slice)

        Returns:
            str: slice of sequence
        """
        return self.sequence[index]

    def __str__(self) -> str:
        """

        Returns:
            str: sequence as a string.
        """
        return self.sequence

    def __repr__(self) -> str:
        """

        Returns:
            str: string that can be used to recreate the object
        """
        return f"{self.__class__.__name__}(sequence='{self.sequence}')"

    @abstractmethod
    def _check_alphabet(self) -> bool:
        """Checks if the sequence contains only allowed characters.

        This method must be implemented by subclasses to validate the sequence
        based on the allowed characters for the specific type of biological
        sequence (e.g., DNA, RNA, protein).

        Returns:
            bool: True if the sequence is valid, False otherwise.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """A base class representing a nucleic acid sequence (DNA or RNA).

    This class provides common functionality for nucleic acid sequences, such as
    complement, reverse, and reverse complement operations. It is intended to be
    subclassed by `DNASequence` and `RNASequence`.

    Attributes:
        _complement_map (dict): A dictionary mapping each nucleotide to its complement.
                                This should be defined in subclasses (e.g., DNA or RNA).


    Raises:
        NotImplementedError: if user tries to create DNA\RNA sequence or find its complement, 
        reverse or reverse_complement through `NucleicAcidSequence` instead of `DNASequence` and `RNASequence`.
    """
    _complement_map = {}

    def _check_alphabet(self) -> bool:
        """Check if the sequence contains only valid nucleotides.

        This method ensures that the sequence contains only characters defined in
        `_complement_map`. It raises a `NotImplementedError` if called directly on
        `NucleicAcidSequence` (use `DNASequence` or `RNASequence` instead).

        Returns:
            bool: True if the sequence is valid, False otherwise.

        Raises:
            NotImplementedError: If called directly on `NucleicAcidSequence`.

        """
        
        if self.__class__ == NucleicAcidSequence:
            raise NotImplementedError(
                "for DNA/RNA creation use DNASequence or RNASequence"
            )
        return set(self.sequence).issubset(self._complement_map.keys())

    def complement(self):
        """This method creates a new sequence where each nucleotide is replaced by its
        complement (as defined in `_complement_map`). It raises a `NotImplementedError`
        if called directly on `NucleicAcidSequence`.

        Raises:
            NotImplementedError: If called directly on `NucleicAcidSequence`.

        Returns:
            NucleicAcidSequence: A new sequence representing the complement
        """
        if self.__class__ == NucleicAcidSequence:
            raise NotImplementedError("")
        return self.__class__(
            "".join(self._complement_map[n] for n in self.sequence)
        )

    def reverse(self):
        """This method creates a new sequence where each nucleotide is replaced by its
        reverse (as defined in `_complement_map`). It raises a `NotImplementedError`
        if called directly on `NucleicAcidSequence`.

        Raises:
            NotImplementedError: If called directly on `NucleicAcidSequence`.

        Returns:
            NucleicAcidSequence: A new sequence representing the reverse of input sequence.
        """
        if self.__class__ == NucleicAcidSequence:
            raise NotImplementedError("")
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        """This method creates a new sequence where each nucleotide is replaced by its
        reverse_complement (as defined in `_complement_map`). It raises a `NotImplementedError`
        if called directly on `NucleicAcidSequence`.

        Raises:
            NotImplementedError: If called directly on `NucleicAcidSequence`.

        Returns:
            NucleicAcidSequence: A new sequence representing the reverse_complement of input sequence.
        """
        if self.__class__ == NucleicAcidSequence:
            raise NotImplementedError("")
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """A class representing a DNA sequence.

    This class inherits from `NucleicAcidSequence` and provides additional functionality
    specific to DNA sequences, such as transcription and complementarity rule.

    Args:
        NucleicAcidSequence (str): DNA string, case sensitive

    Attributes:
        _complement_map (dict): A dictionary mapping each DNA nucleotide to its complement.
                                Includes both uppercase and lowercase mappings.

    """
    _complement_map = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }

    def transcribe(self):
        """method for transcription

        Returns:
            str: transcribed RNA chain as RNASequence
        """
        translation_table = str.maketrans({"T": "U", "t": "u"})
        transcribed_sequence = self.sequence.translate(translation_table)
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
     """A class representing a RNA sequence.

    This class inherits from `NucleicAcidSequence` and provides additional functionality
    specific to DNA sequences, such as complementarity rule.

    Args:
        NucleicAcidSequence (str): RNA string, case sensitive

    Attributes:
        _complement_map (dict): A dictionary mapping each RNA nucleotide to its complement.
                                Includes both uppercase and lowercase mappings.

    """
    _complement_map = {
        "A": "U",
        "G": "C",
        "U": "A",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }


class AminoAcidSequence(BiologicalSequence):
   
    """A class representing an amino acid sequence (protein).
    This class inherits from `BiologicalSequence` and provides functionality
    specific to protein sequences, such as validation and alphabet checking.

    Attributes:
        _valid_amino_acids (set): A set of valid amino acid characters, including
                                  the stop codon symbol as `*`.
        _start_codon (str): The start codon for protein sequences (methionine, "M"). 

    Raises:
        ValueError: in case sequence is empty 
        ValueError: in case characters in sequence are not allowed
        ValueError: in case sequence does not start from methionine 
        ValueError: in case sequence us too short to form a protein

    """
    _valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*")
    _start_codon = "M"

    def _check_alphabet(self) -> bool:
        """Checks if the sequence contains only allowed characters.

        This method must be implemented by subclasses to validate the sequence
        based on the allowed characters for protein).

        Returns:
            bool: True if the sequence is valid, False otherwise.
        
        """
        
        if not self.sequence:
            return False
        return set(self.sequence.upper()).issubset(self._valid_amino_acids)

    def validate_protein(self) -> bool:
        """
        Raises:
            ValueError: in case sequence is empty 
            ValueError: in case characters in sequence are not allowed
            ValueError: in case sequence does not start from methionine 
            ValueError: in case sequence us too short to form a protein

        Returns:
            bool: True if the sequence is a valid protein, otherwise raises an exception
        """
         
        if not self.sequence:
            raise ValueError("Sequence is empty")
        if not self._check_alphabet():
            raise ValueError("")
        if self.sequence[0] != self._start_codon:
            raise ValueError("Sequence does not start with methionine (M)")
        if len(self.sequence) < 20:
            raise ValueError("Sequence is too short to be a valid protein")
        return True


def filter_fastq(
    input_file: str,
    output_file: str = "output_filter.fastq",
    gc_bounds: tuple[float, float] = (0, 100),
    length_bounds: tuple[int, int] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """Filter FASTQ sequences according to length, GC content, and quality thresholds.
    The result is saved in an output file in a directory named 'filtered'.

    Args:
        input_file (str): Path to the input FASTQ file.
        output_file (str, optional): Output file name with extension.
            Defaults to 'output_filter.fastq'.
        gc_bounds (tuple[float, float], optional): GC content bounds (min, max).
            Defaults to (0, 100).
        length_bounds (tuple[int, int], optional): Length bounds (min, max).
            Defaults to (0, 2**32).
        quality_threshold (float, optional): Minimum average quality score.
            Defaults to 0.
    """
    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    filtered_records = []
    for record in SeqIO.parse(input_file, "fastq"):
        if not record.seq:
            continue

        gc_content = gc_fraction(record.seq) * 100
        avg_quality = sum(record.letter_annotations["phred_quality"]) / len(
            record
        )

        if (
            gc_bounds[0] <= gc_content <= gc_bounds[1]
            and length_bounds[0] <= len(record) <= length_bounds[1]
            and avg_quality >= quality_threshold
        ):
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_path, "fastq")
