#Get the COG category of each gene in the model
from IPython import embed
from subprocess import call
import os
from glob import glob
import pandas as pd


def main():
    #Get the core genomes
    core_files = glob('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/*-pan-prots.faa')
    genomes = [element.split('/')[7] for element in core_files]
    from Bio import SeqIO
    #Blast against COGs database    
    ref = '/home/jean-christophe/Documents/Maitrise_UCSD/Comparative_genomics/COG/Blast/prot2003_2014.fa'
    print('gonna run the generate orthologs')
    gen_orthologs(ref, genomes)
    
    #Use blast results to calculate BBHs
    for i in core_files:
        filename = mypath+i
        #out=calc_bbh(ref, i)
        
    #Generate orthogonality matrix
    #bbh_matrix = gen_matrix(ref, core_files)
    embed()
    
def gen_orthologs(ref, genomes): # calculates orthologs for all genomes to a reference (all AAs)
    os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/')
    for g in genomes:
        # genome vs reference
        os.system('makeblastdb -in %s -dbtype prot -out reference'%(ref,))
        outfile= g+'_vs_COG_blast_results'
        print('Now blasting %s against reference, output file: %s'%(g,outfile))
        os.system('blastp -query %s -db reference -outfmt 6 -out %s'%(g,outfile))
        
        # reference vs genome
        #do not compute reference vs genome, irrelevant for that purpose
        '''
        os.system('makeblastdb -in %s -dbtype prot -out genome'%(g,))
        outfile=ref+'_vs_'+g+'_blast_results'
        os.system('blastp -query %s -db genome -outfmt 6 -out %s'%(ref,outfile))
        '''
        
def calc_bbh(query, subject): # calculates BBHs from blast results

    print 'parsing BBHs for', query, subject
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    infile = '%s_vs_%s_blast_results'%(query, subject)
    if not os.path.exists(infile):
        return 
    bbh=pd.read_csv(infile, sep='\t', names=cols)
    bbh2=pd.read_csv('%s_vs_%s_blast_results'%(subject, query), sep='\t', names=cols)
    out = pd.DataFrame()
    for g in bbh.gene.unique():
        res = bbh[bbh.gene==g]
        if len(res)==0:
            continue
        best_hit = res.ix[res.PID.idxmax()].copy()
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene==best_gene]
        if len(res2)==0:
            continue
        best_hit2 = res2.ix[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        # embed()
        if g==best_gene2:
            # print g, '<=>', best_gene
            best_hit['BBH'] = '<=>'
        else:
            # print g, '->', best_gene
            best_hit['BBH'] = '->'
        out=pd.concat([out, pd.DataFrame(best_hit).transpose()])
    out_file = '%s_vs_%s_bbh.csv'%(query, subject)
    out.to_csv(out_file)
    return out  
    
def gen_matrix(ref, genomes): # calculates BBH matrix
    
    #Create an empty dataframe to fill for output
    out = pd.DataFrame()
    
    #For each genome to analyze
    for g in genomes:
        f = '%s_vs_%s_bbh.csv'%(ref, g)
        if not os.path.exists(f):
            continue
        data = pd.read_csv(f, index_col=0)
        data = data[(data.PID>60)&(data.BBH=='<=>')]
        data.index = data.gene
    
        data2 = data[['subject','PID']]
        if len(out)==0:
            out = data2
            out = out.rename(columns={'subject':g,'PID':g+'_PID'})
        else:
            out = pd.merge(out, data2, left_index=True, right_index=True, how='outer')
            out = out.rename(columns={'subject':g,'PID':g+'_PID'})
    return out


if __name__=='__main__':
    main()
