from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def writeFasta(df,filename):
    #Parse the file and make a workable fasta
    my_records = []
    for i,row in df.iterrows():
        #Obtain id, description and sequence from original file
        if i % 2 == 0:
            misc_stuff = str(row.stuff)
            list_of_sutff = misc_stuff.split('|')
            try:
                desc = list_of_sutff[3]
                new_id = list_of_sutff[6].split('^', 1)[0]
            except:
                print('Index out of range for %s file at line %s'%(filename,i))
        elif i % 2 != 0:
            sequence = str(row.stuff)
            #Make SeqRecord element from this
            rec = SeqRecord(Seq(sequence,generic_protein),
                            id =new_id,
                           description = desc)
            my_records.append(rec)

    #Write a proper fasta file from this collection of seqrecord
    from Bio import SeqIO
    SeqIO.write(my_records,filename, 'fasta')
    
'''-------------------------------'''

import pandas as pd
path
os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/')

#Get the filenames in the directory
from os import listdir
from os.path import isfile, join
mypath = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/'
original_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

#Generate new filenames for the fasta
for original_filename in original_files: 
    df = pd.read_table(original_filename,header=-1,names=['stuff'])
    new_filename = original_filename.replace('.txt','.faa')
    writeFasta(df,new_filename)
