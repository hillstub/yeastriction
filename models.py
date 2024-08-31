import os
import subprocess
import tempfile
import json
from typing import List, Optional, Dict, Any, Tuple
from datetime import datetime
from sqlmodel import Field, Relationship, SQLModel, Session, create_engine, select
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import joinedload
from sqlalchemy import JSON, Column
import primer3
import re
import RNA
from dotenv import load_dotenv
from utils import reverse_complement, getRNACentroidStructure

load_dotenv()

DATABASE_URL = os.getenv('DATABASE_URL')

engine = create_engine(DATABASE_URL)
LocalSession = sessionmaker(bind=engine)

MINIMUM_FLANKING_SEQUENCE = 85
KNOCKOUT_LOCUS_PADDING = 60
MINIMUM_KNOCKOUT_LOCUS_LENGTH = 50

class CRISPRSystem(SQLModel, table=True):
    """
    Represents a CRISPR system in the database.

    Attributes:
        id (int): The unique identifier for the CRISPR system.
        name (str): The name of the CRISPR system.
        recognition_sequence_regexp (str): The regular expression for the recognition sequence.
        rna_template_sequence (str): The RNA template sequence for the CRISPR system.
        target_filter_function (str): The name of the function used to filter targets.
        oligo_build_methods (List[Dict]): A list of dictionaries describing oligo build methods.
    """
    id: int = Field(default=None, primary_key=True)
    name: str
    target_filter_function: str
    description: Optional[str]
    recognition_sequence_regexp: str
    rna_template_sequence: str
    oligo_build_methods: Optional[List[Dict[str, Any]]] = Field(sa_column=Column(JSON))


class RNAFold(SQLModel, table=True):
    """
    Represents an RNA fold in the database.

    Attributes:
        id (int): The unique identifier for the RNA fold.
        notation (str): The notation representing the RNA fold structure.
        score (float): The score associated with the RNA fold.
        notation_binding_only (str): The notation representing only the binding part of the RNA fold.
    """
    id: int = Field(default=None, primary_key=True)
    notation: Optional[str]
    notation_binding_only: Optional[str]
    score: Optional[float]

class DiagnosticPrimer(SQLModel, table=True):
    """
    Represents a diagnostic primer in the database.

    Attributes:
        id (int): The unique identifier for the diagnostic primer.
        sequence (str): The sequence of the diagnostic primer.
    """
    id: int = Field(default=None, primary_key=True)
    sequence: Optional[str]

class Strain(SQLModel, table=True):
    """
    Represents a strain in the database.

    Attributes:
        id (int): The unique identifier for the strain.
        name (str): The name of the strain.
        description (str): The description of the strain.
        loci (relationship): Relationship to Locus objects associated with this strain.
    """
    id: int = Field(default=None, primary_key=True)
    name: str = Field(default="", alias='name')
    description: str = Field(default="", alias='description')
    loci: Optional[List["Locus"]] = Relationship(back_populates="strain")

class Target(SQLModel, table=True):
    """
    Represents a target in the database.

    Attributes:
        id (int): The unique identifier for the target.
        sequence (str): The sequence of the target.
        sequence_wo_pam (str): The sequence without the PAM (Protospacer Adjacent Motif).
        position (int): The position of the target in the locus sequence.
        GC_content (float): The GC content of the target sequence.
        z_score (float): The calculated z-score for the target.
        locus_id (int): The ID of the associated locus.
        locus (relationship): Relationship to the associated Locus object.
        crispr_system_id (int): The ID of the associated CRISPR system.
        crispr_system (relationship): Relationship to the associated CRISPRSystem object.
        rna_fold (relationship): Relationship to the associated RNAFold object.
    """
    id: int = Field(default=None, primary_key=True)
    rna_fold_id: int = Field(default=None, foreign_key="rnafold.id")
    rna_fold: Optional["RNAFold"] = Relationship()
    position: Optional[int]
    GC_content: Optional[float]
    sequence_wo_pam: Optional[str]
    sequence: Optional[str]
    locus_id: int = Field(default=None, foreign_key="locus.id")
    locus: Optional["Locus"] = Relationship(back_populates="targets")
    crispr_system_id: int = Field(default=None, foreign_key="crisprsystem.id")
    crispr_system: "CRISPRSystem" = Relationship()
    z_score: Optional[float]

    def get_build_oligos(self, session: Session, dna_build_method: str) -> List[Dict[str, str]]:
        """
        Generate build oligos for the target based on the specified DNA build method.

        Args:
            session (Session): The database session.
            dna_build_method (str): The name of the DNA build method to use.

        Returns:
            List[Dict[str, str]]: A list of dictionaries containing primer names and sequences.
        """
        restricted_env = {
            'target': self.__dict__,
            'reverse_complement': reverse_complement,
            '__builtins__': {}
        }
     
        build_oligos = []

        oligo_build_methods = self.crispr_system.oligo_build_methods


        for oligo_build_method in oligo_build_methods:
            for instruction_build_oligo in oligo_build_method['oligos']:
                if oligo_build_method['name'] == dna_build_method:

                    oligo = {
                        'primer_name': self.locus.display_name + instruction_build_oligo['suffix'],
                        'primer_sequence': eval(instruction_build_oligo['function'], {}, restricted_env)
                    }
                    build_oligos.append(oligo)

        return build_oligos
    
class Locus(SQLModel, table=True):
    """
    Represents a locus in the database.

    Attributes:
        id (int): The unique identifier for the locus.
        orf (str): The Open Reading Frame (ORF) of the locus.
        symbol (str): The symbol associated with the locus.
        sequence (str): The DNA sequence of the locus.
        start_orf (int): The start position of the ORF in the sequence.
        end_orf (int): The end position of the ORF in the sequence.
        strain_id (int): The ID of the associated strain.
        strain (relationship): Relationship to the associated Strain object.
        targets (relationship): Relationship to Target objects associated with this locus.
    """
    id: int = Field(default=None, primary_key=True)
    created: datetime = Field(default_factory=datetime.now)
    sgd_id: str = Field(default='', alias='sgdId')
    orf: str
    symbol: Optional[str]
    strain_id: int = Field(default=None, foreign_key="strain.id")
    strain: Optional["Strain"] = Relationship(back_populates="loci")
    
    sequence: str
    start_orf: int = Field(..., alias='startOrf')
    end_orf: int = Field(..., alias='endOrf')
    
    forward_diagnostic_primer_id: Optional[int] = Field(default=None, foreign_key="diagnosticprimer.id")
    forward_diagnostic_primer: Optional["DiagnosticPrimer"] = Relationship(sa_relationship_kwargs={"foreign_keys": "Locus.forward_diagnostic_primer_id"})
    
    reverse_diagnostic_primer_id: Optional[int] = Field(default=None, foreign_key="diagnosticprimer.id")
    reverse_diagnostic_primer: Optional["DiagnosticPrimer"] = Relationship(sa_relationship_kwargs={"foreign_keys": "Locus.reverse_diagnostic_primer_id"})
    
    targets: List["Target"] = Relationship(back_populates="locus")

    @property
    def display_name(self) -> str:
        return self.symbol if self.symbol else self.orf

    @property
    def repair_oligo_fw(self) -> str:
        return self.sequence[self.start_orf - 60:self.start_orf] + self.sequence[self.end_orf:self.end_orf + 60]

    @property
    def repair_oligo_rv(self) -> str:
        return reverse_complement(self.repair_oligo_fw)

    def get_diagnostic_primers(self, session: Session) -> Tuple[str, str]:
        """
        Generate diagnostic primers for the locus.

        Args:
            session (Session): The database session.

        Returns:
            Tuple[str, str]: A tuple containing the forward and reverse diagnostic primer sequences.
        """
        def is_valid_dna_sequence(sequence: str) -> bool:
            """
            Check if the given sequence is a valid DNA sequence.

            Args:
                sequence (str): The DNA sequence to check.

            Returns:
                bool: True if the sequence is valid, False otherwise.
            """
            return bool(re.match("^[AGTC]*$", sequence))
                
        if self.forward_diagnostic_primer_id and self.reverse_diagnostic_primer_id:
            return self.forward_diagnostic_primer.sequence, self.reverse_diagnostic_primer.sequence

        if self.start_orf <= MINIMUM_FLANKING_SEQUENCE or (len(self.sequence) - self.end_orf) <= MINIMUM_FLANKING_SEQUENCE:
            return '', ''

        ko_locus = self.sequence[:self.start_orf-KNOCKOUT_LOCUS_PADDING] + self.sequence[self.end_orf+KNOCKOUT_LOCUS_PADDING:]
        if len(ko_locus) <= MINIMUM_KNOCKOUT_LOCUS_LENGTH:
            return '', ''
        
        # we don't support ambigious bases.
        if not re.match("^[AGTC]*$", ko_locus):
            return '', ''

        input_sequence = ko_locus
        seq_args = {
            'SEQUENCE_TEMPLATE': input_sequence,
        }
        global_args = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_NUM_RETURN': 1,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [
                [0, self.start_orf - 60,
                self.start_orf + 60, len(input_sequence) - (self.start_orf + 60)]
            ],
            'PRIMER_PRODUCT_SIZE_RANGE': [[250, 750]],
            'PRIMER_GC_CLAMP': 1
        }

        primer3_results = primer3.bindings.design_primers(seq_args, global_args)

        forward_primer_sequence = primer3_results.get('PRIMER_LEFT_0_SEQUENCE', '')
        reverse_primer_sequence = primer3_results.get('PRIMER_RIGHT_0_SEQUENCE', '')

        if forward_primer_sequence and reverse_primer_sequence:
            forward_primer = DiagnosticPrimer(sequence=forward_primer_sequence)
            reverse_primer = DiagnosticPrimer(sequence=reverse_primer_sequence)
            session.add(forward_primer)
            session.add(reverse_primer)
            session.commit()

            self.forward_diagnostic_primer_id = forward_primer.id
            self.reverse_diagnostic_primer_id = reverse_primer.id
            session.commit()

        return forward_primer_sequence, reverse_primer_sequence
    

def get_locus_from_database(session: Session, locus_id: int, crispr_system_id: int) -> Locus:
    """
    Get the locus from the database with the given ID and CRISPR system ID.

    Args:
        session (Session): The database session.
        locus_id (int): The ID of the locus.
        crispr_system_id (int): The ID of the CRISPR system.

    Returns:
        Locus: The locus object with the targets and RNA fold loaded.
    """
    # Query the locus with targets and rna_fold loaded
    locus = session.query(Locus).options(
        joinedload(Locus.targets.and_(Target.crispr_system_id == crispr_system_id)).joinedload(Target.rna_fold)
    ).filter(Locus.id == locus_id).first()

    return locus

def no_filter(session: Session, locus: Locus, targets: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    A no-op filter function that returns the input targets without modification.

    Args:
        session (Session): The database session (unused).
        locus (Locus): The Locus object containing the targets (unused).
        targets (List[Dict[str, Any]]): A list of target dictionaries.

    Returns:
        List[Dict[str, Any]]: The input list of target dictionaries, unmodified.
    """
    return targets

def filter_cas9_targets_with_bowtie(session: Session, locus: Locus, targets: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Filter Cas9 targets using Bowtie alignment.

    Args:
        session (Session): The database session.
        locus (Locus): The Locus object containing the targets.
        targets (List[Dict[str, Any]]): A list of target dictionaries to filter.

    Returns:
        List[Dict[str, Any]]: A filtered list of target dictionaries.
    """
    # Add 8 different sequences per target
    nucleotides = ['A', 'T', 'G', 'C']
    variants = []

    for index, target in enumerate(targets):
        for variant_index, nucleotide in enumerate(nucleotides):
            variants.append(f">target_{index}_{variant_index * 2}\n{target['sequence_wo_pam'] + nucleotide + 'GG'}")
            variants.append(f">target_{index}_{variant_index * 2 + 1}\n{target['sequence_wo_pam'] + nucleotide + 'AG'}")

    # Create a temporary file to store the variants
    with tempfile.NamedTemporaryFile(mode='w+') as tmpfile:
        tmpfile.write('\n'.join(variants))
        tmpfile_path = tmpfile.name

        # Call the bowtie executable
        genome_path = os.path.join(os.getenv('GENOMES_DIR'), locus.strain.name)
        genome_path = locus.strain.name
        bowtie_command = f'bowtie -k 2 -v 3 {genome_path} --suppress 2,3,4,5,6,7,8 -f {tmpfile_path} 2> /dev/null | uniq -c | awk \'{{print $2,$1}}\''
        result = subprocess.run(bowtie_command, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Bowtie error: {result.stderr}")
            return []

        bowtie_hits = result.stdout.split('\n')

        indices_with_multiple_hits = []
        for hit in bowtie_hits:
            if hit:
                row = hit.split(' ')
                variant_id = row[0].split('_')
                target_index = int(variant_id[1])
                if int(row[1]) >= 2:
                    indices_with_multiple_hits.append(target_index)

        # Filter the targets based on bowtie results
        filtered_targets = [targets[index] for index in range(len(targets)) if index not in indices_with_multiple_hits]

        return filtered_targets


def search_targets(session: Session, locus_id: int, crispr_system_id: int) -> List[Target]:
    """
    Search for targets in a given locus for a specific CRISPR system.

    Args:
        session (Session): The database session.
        locus_id (int): The ID of the locus to search in.
        crispr_system_id (int): The ID of the CRISPR system to use.

    Returns:
        List[Target]: A list of Target objects found for the given locus and CRISPR system.
    """
    crispr_system = session.query(CRISPRSystem).filter(CRISPRSystem.id == crispr_system_id).first()
    
    locus = get_locus_from_database(session=session, locus_id=locus_id, crispr_system_id=crispr_system_id)

    if not locus:
        return []

    if locus.targets:
        return locus.targets

    reg = re.compile(rf'{crispr_system.recognition_sequence_regexp}')
    orf_sequence = locus.sequence[locus.start_orf:locus.end_orf]
    targets = []

    for match in reg.finditer(orf_sequence):
        targets.append({
            'sequence': match.group('target_sequence_with_pam'),
            'sequence_wo_pam': match.group('target_sequence_without_pam')
        })
    
    for match in reg.finditer(reverse_complement(orf_sequence)):
        targets.append({
            'sequence': match.group('target_sequence_with_pam'),
            'sequence_wo_pam': match.group('target_sequence_without_pam')
        })

    targets = [target for target in targets if 'TTTTTT' not in target['sequence_wo_pam']]
    
    rna_template_sequence = crispr_system.rna_template_sequence
    position_target_sequence_in_template_sequence = rna_template_sequence.format(target_sequence_without_pam='*').find('*')
    target_objects = []
    gc_contents = []
    rna_fold_scores = []

    if crispr_system.target_filter_function == 'filter_cas9_targets_with_bowtie':
        targets = filter_cas9_targets_with_bowtie(session=session, locus=locus, targets=targets)
    elif crispr_system.target_filter_function == 'no_filter':
        targets = no_filter(session=session, locus=locus, targets=targets)
    else:
        raise ValueError('No filtering function defined!')

    for target in targets:
        target['position'] = orf_sequence.find(target['sequence']) if target['sequence'] in orf_sequence else orf_sequence.find(reverse_complement(target['sequence']))
        target['GC_content'] = (target['sequence_wo_pam'].count('G') + target['sequence_wo_pam'].count('C')) / len(target['sequence_wo_pam'])
        rna_sequence = rna_template_sequence.format(target_sequence_without_pam=target['sequence_wo_pam'])
        structure, score = getRNACentroidStructure(rna_sequence)
        notation_binding_only = structure[position_target_sequence_in_template_sequence:position_target_sequence_in_template_sequence+len(target['sequence_wo_pam'])]
        score = notation_binding_only.count('.') / len(target['sequence_wo_pam'])
        target['rna_fold'] = {'notation': structure, 'score': score, 'notation_binding_only': notation_binding_only}

        gc_contents.append(target['GC_content'])
        rna_fold_scores.append(score)

        # Calculate ranges
    ranges = {
        'rna_fold_score': {'min': min(rna_fold_scores), 'max': max(rna_fold_scores)},
        'gc_content': {'min': min(gc_contents), 'max': max(gc_contents)}
    }

    for target in targets:
        target_score = 0
        target_score += 1 - ((target['GC_content'] - ranges['gc_content']['min']) / (ranges['gc_content']['max'] - ranges['gc_content']['min']))
        target_score += (target['rna_fold']['score'] - ranges['rna_fold_score']['min']) / (ranges['rna_fold_score']['max'] - ranges['rna_fold_score']['min'])
        
        target['z_score'] = target_score

        new_rna_fold = RNAFold(notation=target['rna_fold']['notation'], score=target['rna_fold']['score'], notation_binding_only=target['rna_fold']['notation_binding_only'])
        new_target = Target(sequence=target['sequence'], sequence_wo_pam=target['sequence_wo_pam'],
                            position=target['position'], GC_content=target['GC_content'], z_score = target['z_score'],
                            locus_id=locus.id, crispr_system_id=crispr_system.id, rna_fold=new_rna_fold)
        target_objects.append(new_target)
        session.add(new_rna_fold)
        session.add(new_target)
    session.commit()

    locus = get_locus_from_database(session=session, locus_id=locus_id, crispr_system_id=crispr_system_id)
    return locus.targets

def initialize_database():
    """
    Initialize the database by creating the tables and inserting initial data.
    """
    # Create the tables if the database doesn't exist
    if not os.path.exists(DATABASE_URL.replace('sqlite:///','')):
        SQLModel.metadata.create_all(engine)

        # Insert initial data into crisprsystem

        initial_data = CRISPRSystem(
            id=1,
            name="Cas9 (NGG)",
            description="CRISPR system with S. pyogenes Cas9 nuclease",
            recognition_sequence_regexp="(?P<target_sequence_with_pam>(?P<target_sequence_without_pam>[ATGC]{20})(?P<pam_sequence>[ATGC]GG))",
            oligo_build_methods=[{"name": "pROS", "oligos": [{"suffix": " pROS fw", "function": "\"tgcgcatgtttcggcgttcgaaacttctccgcagtgaaagataaatgatc\"+target[\"sequence_wo_pam\"]+\"gttttagagctagaaatagcaagttaaaataag\""}]}, {"name": "pMEL", "oligos": [{"suffix": " pMEL fw", "function": "\"tgcgcatgtttcggcgttcgaaacttctccgcagtgaaagataaatgatc\"+target[\"sequence_wo_pam\"]+\"gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaac\""}, {"suffix": " pMEL rv", "function": "reverse_complement(\"tgcgcatgtttcggcgttcgaaacttctccgcagtgaaagataaatgatc\"+target[\"sequence_wo_pam\"]+\"gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaac\")"}]}],
            rna_template_sequence='{target_sequence_without_pam}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTT',
            target_filter_function='filter_cas9_targets_with_bowtie'
        )

        with LocalSession() as session:
            session.add(initial_data)
            session.commit()