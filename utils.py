import RNA

def reverse_complement(sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        sequence (str): The input DNA sequence.

    Returns:
        str: The reverse complement of the input sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def getRNACentroidStructure(sequence: str, temperature: float = 30.) -> tuple[str, float]:
    """
    Calculate the centroid structure of an RNA sequence.

    Args:
        sequence (str): The input RNA sequence.
        temperature (float): The temperature at which the RNA structure is calculated. Default is 30.0.
    Returns:
        tuple[str, float]: A tuple containing the centroid structure (str) and its score (float).
    """
    settings = RNA.md()
    settings.temperature = temperature
    settings.dangles = 2
    settings.noLP = 1

    fc_obj = RNA.fold_compound(sequence, settings)
    fc_obj.pf()
    structure, value = fc_obj.centroid()
    return structure, value
