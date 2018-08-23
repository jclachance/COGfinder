import pandas as pd
import os

#COG conversion files (constants)
CONV_COG = pd.read_csv('/home/jean-christophe/Documents/Maitrise_UCSD/Comparative_genomics/COG/Blast/cog2003-2014.csv',
                       names=['gi','specie','stuff','stuff2','stuff3', 'COG','binary'])
LETTERS_COG = pd.read_table('/home/jean-christophe/Documents/Maitrise_UCSD/Comparative_genomics/COG/Blast/cognames2003-2014.tab')
UNIQ_LETTERS = pd.read_table('/home/jean-christophe/Documents/Maitrise_UCSD/Comparative_genomics/COG/Blast/fun2003-2014.tab')

#Sort data
#Read the dataframe of blast results
def trimBlast(filename):
    TRESHOLD = 1e-10
    os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/Blast_results/Pan')
    df = pd.read_table(filename,
                      names = ['Mflorum_prot','db_prot','coverage','start','end','stuff','stuff2','stuff3',
                              'stuff4','stuff5','e_val','stuff6'])
    df_trim = df[df.e_val<=TRESHOLD]
    return df_trim
'''---------------------------------'''

def extractGI(df):
    Mflorum_prot, gi_db = [], []
    for i,row in df.iterrows():
        Mflorum_prot.append(row.Mflorum_prot)
        prot_db =  row.db_prot.split('|')
        gi = prot_db[1]
        gi_db.append(gi)
    #Make dataframe
    df2 = pd.DataFrame({'Mflorum_prot': Mflorum_prot,
                   'prot_db':gi_db})
    df2.prot_db = pd.to_numeric(df2.prot_db)
    return df2
'''---------------------------------'''

def makeFlorumCOG(gi_df):
    PATH = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/'
    #Merge conversion to COGcategories and blast results
    merged_df = pd.merge(left=gi_df, right=CONV_COG,
                        how='inner',left_on='prot_db',
                        right_on='gi')
    #Get only the columns that matter
    temp_df = pd.DataFrame({'Mflorum_prot':merged_df.Mflorum_prot,
                             'COG':merged_df.COG})
    #Eliminate duplicates
    temp_df1 = temp_df.drop_duplicates()
    #Group by COG categories
    grouped = temp_df1.groupby(temp_df1.COG)
    pre_florumCOG = grouped.agg(lambda x: ' and '.join(x))
    filename = PATH + 'temp.csv'
    pre_florumCOG.to_csv(filename)
    florumCOG = pd.read_csv(filename,names=['COG','Mflorum_prot'],skiprows=2)
    return florumCOG
'''---------------------------------'''
   
def COGletters(florumCOG):
    PATH = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/'
    lettersCOG = pd.merge(left=LETTERS_COG,right= florumCOG,
                         left_on='# COG',right_on='COG')
    grouped = lettersCOG.groupby('func')
    temp_df = grouped.agg(lambda x: ' and '.join(x))
    filename = PATH +'temp.csv'
    temp_df.to_csv(filename)
    temp_df1 = pd.read_csv(filename,
                           names=['letter','func', 'name', 'COG', 'Mflorum_prot'],
                           skiprows=1)
    del temp_df1['func']
    return temp_df1
'''---------------------------------'''

def funcCategories(lettersCOG):
    PATH = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/'
    genes, codes = [],[]
    for i, row in lettersCOG.iterrows():
        for l in row.letter:
            for code in UNIQ_LETTERS['# Code']:
                if l == code:
                    codes.append(code)
                    genes.append(row.Mflorum_prot)
    temp_df = pd.DataFrame({'Functional_Categories':codes,'Mflorum_prot':genes})
    grouped = temp_df.groupby('Functional_Categories')
    temp_df1 = grouped.agg(lambda x:' and '.join(x))
    filename = PATH+'temp.csv'
    temp_df1.to_csv(filename)
    final_df = pd.read_csv(filename,skiprows=1,
                          names= ['Functionnal_Categories','Mflorum_prot'])
    return final_df
'''---------------------------------'''

def writeCSV(filename,letters_and_genes):
    os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/')
    strain = filename.rsplit('-pan',1)[0]
    new_filename = strain + 'pan_COG.csv'
    letters_and_genes.to_csv(new_filename)

def main():
    #Extract blast_results from previous cells
    #Blast results have been moved via unix command in terminal
    #Optionally run this command to move data previosuly generated to the folder to work in
    #os.system('mv *blast_results /home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/Blast_results/')
    os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/Blast_results/Pan')
    #Generate list of Blast results files to work with
    from os import listdir
    from os.path import isfile,join
    mypath = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/Blast_results/Pan'
    blast_results = [f for f in listdir(mypath) if isfile(join(mypath, f)) and 'blast_results'in f]
    #Process the data to associate with COG categories
    
    for filename in blast_results:
        df = trimBlast(filename)
        #Extract gi
        gi_df = extractGI(df)
        #Merge with CONV_COG
        florumCOG = makeFlorumCOG(gi_df)
        #Associate with letters 
        lettersCOG = COGletters(florumCOG)
        #Match proteins with functional categories
        letters_and_genes = funcCategories(lettersCOG)
        #Write to csv file
        writeCSV(filename, letters_and_genes)
        
    print(letters_and_genes) 
    
        
if __name__=='__main__':
    main()
