#Generate Histogram for pan genomes
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pylab import rcParams
from os import listdir
from os.path import isfile,join

#Make graphs
def number_of_genes(df):
    number_of_genes = []
    for i,row in df.iterrows():
        mini_list = row.Mflorum_prot.split(' and ')
        number_of_genes.append(len(mini_list))
    new_df = pd.DataFrame({'Functional_Categories':df.Functionnal_Categories,
                           'Number_of_genes':number_of_genes})
    sorted_new_df = new_df.sort_values('Number_of_genes',ascending = False)
    return sorted_new_df  

def categoryLetter(df):
    PATH = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/'
    df2 = pd.merge(left=UNIQ_LETTERS, right=df, left_on='# Code',
                    right_on='Functional_Categories', how='inner')
    del df2['Functional_Categories']
    df2['Functional_Categories'] = df2[['# Code', 'Name']].apply(lambda x: ' '.join(x),axis=1)
    del df2['# Code']
    del df2['Name']
    df2.sort_values('Number_of_genes',inplace=True)
    filename = PATH + 'temp.csv'
    df2.to_csv(filename)
    df3 = pd.read_csv(filename)
    del df3['Unnamed: 0']
    return df3

def plotHistogram(filename,df):
    PATH = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/Figures/'
    
    plt.clf()
    #Get labels
    y_labels = []
    for r in df.Functional_Categories:
        y_labels.append(r)
    y_pos = np.arange(df.Functional_Categories.size)
    width = 0.5
    #Set figure size
    fig = plt.figure()
    default_size = fig.get_size_inches()
    fig.set_size_inches(15,10)
    #Plot a histogram
    #Set color palette with seaborn
    color_palette = sns.color_palette('GnBu_d',len(y_labels))

    plt.barh(y_pos, df.Number_of_genes, align='center', alpha=0.6,height=1.0,
             linewidth=0.5, color=color_palette)
    #Set figure display parameters
    plt.axis([0,200,-1,20])
    plt.yticks(y_pos,y_labels,fontsize=12)
    plt.ylabel('Functional categories',fontsize=24)
    plt.xlabel('Number of genes',fontsize=24)
    #Save figure
    fig_filename = PATH + filename.rsplit('.csv',1)[0] +'.svg'
    print(fig_filename)
    plt.savefig(fig_filename)
    
def plot():

    mypath = '/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG/'
    csv_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f.endswith('pan_COG.csv')]
    os.chdir('/home/jean-christophe/Documents/Maitrise_UCSD/MultiStrain_meso/Data/COG')
    for filename in csv_files:
        print(filename)
        df = pd.read_csv(filename)
        del df['Unnamed: 0']
        number_per_category = number_of_genes(df)
        #Make dataframe with sorted categories based on number of genes in it
        categoryDF = categoryLetter(number_per_category)
        plotHistogram(filename, categoryDF)


