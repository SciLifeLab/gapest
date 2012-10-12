'''
Created on Sep 23, 2011

@author: ksahlin
'''
#import numpy

#import Contig,Scaffold

def ReadInContigseqs(contigfile):
    cont_dict={}
    k=0
    temp=''
    accession=''
    for line in contigfile:
        if line[0]=='>' and k==0:
            accession=line[1:].strip()
            cont_dict[accession]=''
            k+=1
        elif line[0]=='>':
            cont_dict[accession]=temp
            temp=''
            accession=line[1:].strip()
        else:
            temp+=line.strip()
    cont_dict[accession]=temp
    return(cont_dict)

def Main(contigfile_,bamfile,mean,edge_support,read_len,ratio,std_dev):
    from time import time    
    tot_start=time()
#create list of list f that will print all scaffolds, contiger in 
    F=[] #list of (ordered) lists of tuples containing (contig_name, orientation, gap_size) tuple is contig and list is scaffold

    Scaffolds={}     #scaffold dict with contig objects for easy fetching of all contigs in a scaffold
# global indicator for scaffolds, used to index scaffolds when they are created
    scaffold_indexer=1
# contig dict that stores contig objects
    Contigs={}
    #print cont_theshold
        
    import CreateGraph as CG
    import GapCalculator as GC
    #Read in the sequences of the contigs in memory
    contigfile=open(contigfile_,'r')
    C_dict=ReadInContigseqs(contigfile)
#iterate over libraries

    start=time()
    if read_len:
        read_length=read_len
    else: #if not specified we set the read length to 50 by default
        read_length=50
    if std_dev:
        sigma=std_dev
    else:
        sigma=False
    (G,Contigs,Scaffolds,F,scaffold_indexer)=CG.PE(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_length)      #Create graph, single out too short contigs/scaffolds and store 
            
    GC.GapEstimator(G,Contigs,Scaffolds,mean,sigma,read_length,edge_support)





if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    
    parser = OptionParser()

    parser.add_option("-s", "--stddev", dest="stddev", default=0,
                      help="estimated standard deviation of libraries", type="int")
    
    parser.add_option("-f", "--bamfile", dest="bamfiles",
                      help="Name of bamfile", type="string")
    
    parser.add_option("-r",dest="readlen",
                      help="Mean read length ",type="int")
    
    parser.add_option("-m",dest="mean",
                      help="mean of insert library",type="int")
    
    parser.add_option("-w",dest="relweight", default=3,
                      help="treshold value for the relative weight of an edge",type="int")
    
    
    parser.add_option("-c", dest="contigfile",
                      help="file of contigs",type="string")
    
    parser.add_option("-e",dest="edgesupport", nargs=1, default=10,
                      help="treshold value for the least nr of links that is needed to create an edge ",type="int")
       
    (options, args) = parser.parse_args()       
        


    #options.qacomputefile not needed yet, is not implemented yet
    Main(options.contigfile,options.bamfiles,options.mean,options.edgesupport,options.readlen,options.relweight,options.stddev)
        
        
