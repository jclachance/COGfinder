from IPython import embed
from subprocess import call
import os
from glob import glob
import pandas as pd

def _get_reference():
    import os.path
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    file_path = os.path.join(parent_dir, 'data/cog2003-2014.csv')
    cog_reference = pd.read_csv(file_path, )

    return cog_reference

def _generate_orthologs(ref, genomes, outdir):
    """
    Calculates orthologs for all genomes to a reference

    :param ref: The COG reference
    :param genomes: Genomes on which to find COG categories
    :return:
    """
    os.chdir(outdir)
    for g in genomes:
        # genome vs reference
        os.system('makeblastdb -in %s -dbtype prot -out reference' % (ref,))
        outfile= g+'_vs_COG_blast_results'
        print('Now blasting %s against reference, output file: %s' % (g,outfile))
        os.system('blastp -query %s -db reference -outfmt 6 -out %s' % (g,outfile))


def blast(path_to_fasta,outdir):
    """
    This function takes the path to a folder containing all fasta files to be blasted against the COG database.
    Note that a write_fasta function is available in COGfinder/util.

    :param path_to_fasta:
    :param outdir:
    :return:
    """
    #1- Extract the genomes from the folder containing the fasta files
    core_files = glob(path_to_fasta)
    genomes = [element.split('/')[7] for element in core_files]

    #2- Blast against COGs database
    from Bio import SeqIO
    #Getting the reference COG database (provided with package)

    REFERENCE = _get_reference()
    print('Generating orthologs')
    _generate_orthologs(REFERENCE, genomes, outdir)



