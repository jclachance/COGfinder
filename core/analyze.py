import pandas as pd
import os
from os import listdir
from os.path import isfile,join

def _get_reference():
    import os.path
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    file_path = os.path.join(parent_dir, 'data/cog2003-2014.csv')
    cog_reference = pd.read_csv(file_path,
                                names=['gi', 'specie', 'stuff', 'stuff2', 'stuff3', 'COG', 'binary'])

    return cog_reference

def _get_conversion(filename):
    import os.path
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    file_path = os.path.join(parent_dir, filename)

    return pd.read_table(file_path)

#Sort data
#Read the dataframe of blast results
def _trim_blast(filename,result_path, THRESHOLD):

    os.chdir(result_path)
    df = pd.read_table(filename,
                      names = ['target_prot','db_prot','coverage','start','end','stuff','stuff2','stuff3',
                              'stuff4','stuff5','e_val','stuff6'])
    df_trim = df[df.e_val<=TRESHOLD]

    return df_trim
'''---------------------------------'''

def _extract_gene_identifier(df):
    target_prot, gi_db = [], []
    for i,row in df.iterrows():
        target_prot.append(row.target_prot)
        prot_db = row.db_prot.split('|')
        gi = prot_db[1]
        gi_db.append(gi)
    #Make dataframe
    df2 = pd.DataFrame({'target_prot': target_prot,
                   'prot_db':gi_db})
    df2.prot_db = pd.to_numeric(df2.prot_db)
    return df2


def _get_target_COG(gi_df,outdir,CONV_COG):
    PATH = outdir
    #Merge conversion to COGcategories and blast results
    merged_df = pd.merge(left=gi_df, right=CONV_COG,
                        how='inner',left_on='prot_db',
                        right_on='gi')

    #Get only the columns that matter
    temp_df = pd.DataFrame({'target_prot':merged_df.target_prot,
                             'COG':merged_df.COG})

    #Eliminate duplicates
    temp_df1 = temp_df.drop_duplicates()

    #Group by COG categories
    grouped = temp_df1.groupby(temp_df1.COG)
    pre_target_COG = grouped.agg(lambda x: ' and '.join(x))
    filename = PATH + 'temp.csv'
    pre_target_COG.to_csv(filename)
    target_COG = pd.read_csv(filename,names=['COG','target_prot'],skiprows=2)
    return target_COG

   
def _get_COG_letters(target_COG, outdir, LETTERS_COG):
    PATH = outdir
    lettersCOG = pd.merge(left=LETTERS_COG,right= target_COG,
                         left_on='# COG',right_on='COG')
    grouped = lettersCOG.groupby('func')
    temp_df = grouped.agg(lambda x: ' and '.join(x))
    filename = PATH +'temp.csv'
    temp_df.to_csv(filename)
    temp_df1 = pd.read_csv(filename,
                           names=['letter','func', 'name', 'COG', 'target_prot'],
                           skiprows=1)
    del temp_df1['func']
    return temp_df1


def _get_functional_categories(lettersCOG, outdir, UNIQ_LETTERS):
    PATH = outdir
    genes, codes = [],[]
    for i, row in lettersCOG.iterrows():
        for l in row.letter:
            for code in UNIQ_LETTERS['# Code']:
                if l == code:
                    codes.append(code)
                    genes.append(row.Mflorum_prot)

    temp_df = pd.DataFrame({'Functional_Categories':codes,'target_prot':genes})
    grouped = temp_df.groupby('Functional_Categories')
    temp_df1 = grouped.agg(lambda x:' and '.join(x))
    filename = PATH + 'temp.csv'
    temp_df1.to_csv(filename)
    final_df = pd.read_csv(filename,skiprows=1,
                          names= ['Functionnal_Categories','target_prot'])
    return final_df


def _write_csv(filename,letters_and_genes):
    os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/')
    strain = filename.rsplit('-pan',1)[0]
    new_filename = strain + 'pan_COG.csv'
    letters_and_genes.to_csv(new_filename)


def analyze(result_path, outdir, THRESHOLD=1e-10):
    """
    Blast results from the previous steps are filtered and COG categories associated with each genes are analyzed.

    :param result_path: Path to folder containing blast results
    :param outdir: Path to folder to write to
    :param THRESHOLD: e-value THRESHOLD for protein match
    """
    #1- Get all COG conversion files
    CONV_COG = _get_reference()
    LETTERS_COG = _get_conversion('data/cognames2003-2014.tab')
    UNIQ_LETTERS = _get_conversion('data/fun2003-2014.tab')

    #2- Extract blast_results from blast.py
    # Generate list of Blast results files to work with
    mypath = result_path
    blast_results = [f for f in listdir(mypath) if isfile(join(mypath, f)) and 'blast_results'in f]

    #3- Process the data to associate genes with COG categories
    for filename in blast_results:
        df = _trim_blast(filename,result_path,THRESHOLD)
        #Extract gi
        gi_df = _extract_gene_identifier(df)
        #Merge with CONV_COG
        florumCOG = _get_target_COG(gi_df,outdir,CONV_COG)
        #Associate with letters
        lettersCOG = _get_COG_letters(florumCOG,outdir, LETTERS_COG)
        #Match proteins with functional categories
        letters_and_genes = _get_functional_categories(lettersCOG,outdir,UNIQ_LETTERS)
        #Write to csv file
        _write_csv(filename, letters_and_genes)

    
        

