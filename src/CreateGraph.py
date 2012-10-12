'''
Created on Sep 23, 2011

@author: ksahlin
'''

import pysam
import Contig,Scaffold
import networkx as nx
from collections import defaultdict

def PE(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len):
    G=nx.Graph()
    #print 'Parsing BAM file...'
    #read_len=50
    #informative_pair={81:(False,True),97:(True,False),113:(False,False),65:(True,True)}
    #I switched to look at mates instead since BWA can give false flag combinations for
    # read-mate when read is mapped but not mate eg 97-149 81-165. But the reverse
    #does not happen.
    #informative_pair={161:(True,False),145:(False,True),129:(True,True),177:(False,False)}
    #threshold=800
    def AddEdges(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len):
        #Clean contig_library
        singeled_out=0
        cont_lengths= bam_file.lengths
        cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
        #print cont_lengths
        cont_names = bam_file.references
        ####### WHEN ADDING SHORTER CONTIGS NOT INCLUDED IN THE SCAFFOLDING, 
        ####### WE NEED TO ALSO INITIALIZE OBJECTS FOR THESE, THIS SHOULD BE DONE SOMEWHERE HERE
        for i in range(0,len(cont_names)):
            if cont_lengths[i] >= 300:
                C=Contig.contig(cont_names[i])   # Create object contig
                C.length = cont_lengths[i]
                C.scaf_length = C.length        # Initially, scaffold consists of only this contig
                C.direction = True              # always in same direction first, False=reverse
                C.position = 0                  #position always 0
                C.links = {}
                Contigs[C.name] = C              # Create a dict with name as key and the object container as value
                S=Scaffold.scaffold('s'+str(scaffold_indexer),[C],C.length)  # Create object scaffold
                Scaffolds[S.name]=S
                C.scaffold=S.name
                G.add_node((S.name,'L'),length=cont_lengths[i])
                G.add_node((S.name,'R'),length=cont_lengths[i])
                scaffold_indexer+=1
        
        
        #Create "node graph" of contigs (that passed the length criteria). Each having a left and right node
        #print 'Nr of contigs/scaffolds included in scaffolding: '+ str(len(Scaffolds))#,Scaffolds.keys()
        
        for scaffold_ in Scaffolds:
            G.add_edge((scaffold_,'L'),(scaffold_,'R'),nr_links=None)    #this is a scaffold object but can be both a single contig or a scaffold.
        
        
        # Create the link edges in the graph by fetching info from bam file
        
        fishy_edges = defaultdict(int)
        for alignedread in bam_file:
            try: #check that read is aligned OBS: not with is_unmapped since this flag is fishy for e.g. BWA
                contig1=bam_file.getrname(alignedread.rname)
                contig2=bam_file.getrname(alignedread.mrnm)
            except ValueError:
                continue  
            if contig1 in Contigs and contig2 in Contigs:
                #TODO: this if-statement is an ad hoc implementation to deal with BWA's buggy SAM-flag reporting
                #if BWA fixes this -> remove this statement. If the links in fishy edges is equal to or ore than
                #the links in the graph G or G'. The edge will be removed.
                if alignedread.is_unmapped and alignedread.is_read1: # and contig1 != contig2: 
                    #Some BWA error in mappings can still slip through, these edges are caracterized by very few links                 
                    cont_obj1 = Contigs[contig1]
                    scaf_obj1 = Scaffolds[cont_obj1.scaffold]
                    cont_obj2 = Contigs[contig2]
                    scaf_obj2 = Scaffolds[cont_obj2.scaffold]
                    
                    if scaf_obj2.name != scaf_obj1.name:
                        (side1,side2) = CheckDir(cont_obj1,cont_obj2,alignedread) 
                        #get scaffold name for contig
                        s1 = Contigs[contig1].scaffold #if contig1 in Contigs else small_contigs[contig1].scaffold
                        s2 = Contigs[contig2].scaffold #if contig2 in Contigs else small_contigs[contig2].scaffold   
                        fishy_edges[((s1,side1),(s2,side2))] +=1
                        fishy_edges[((s2,side2),(s1,side1))] +=1
                
                #if contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
                if contig1 != contig2 and alignedread.is_read2 and not alignedread.is_unmapped and alignedread.mapq  > 20:
                    (read_dir,mate_dir) = (not alignedread.is_reverse,not alignedread.mate_is_reverse )
                    scaf1=Contigs[contig1].scaffold
                    scaf2=Contigs[contig2].scaffold                    
                    #Calculate actual position on scaffold here
                    #position1 cont/scaf1
                    cont_dir1 = Contigs[contig1].direction  #if pos : L if neg: R
                    cont1_pos = Contigs[contig1].position
                    readpos = alignedread.pos
                    cont1_len = Contigs[contig1].length
                    s1len = Scaffolds[scaf1].s_length
                    #position1 cont1/scaf1                        
                    cont_dir2 = Contigs[contig2].direction
                    cont2_pos = Contigs[contig2].position
                    matepos = alignedread.mpos
                    cont2_len = Contigs[contig2].length
                    s2len = Scaffolds[scaf2].s_length 
                    (obs,scaf_side1,scaf_side2)=PosDirCalculatorPE(cont_dir1,read_dir,cont1_pos,readpos,s1len,cont1_len,cont_dir2,mate_dir,cont2_pos,matepos,s2len,cont2_len,read_len) 
                    if obs < mean+ 6*std_dev: 
                        if (scaf2,scaf_side2) not in G[(scaf1,scaf_side1)]:
                            G.add_edge((scaf2,scaf_side2),(scaf1,scaf_side1),nr_links=1,gap_dist=[obs])
                        #print 'Added edge'
                        else:
                            G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['nr_links'] += 1
                            #print 'edge'
                            G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['gap_dist'].append(obs)                         
            
            elif contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
                ########################Use to validate scaffold herein previous step here
                pass
        RemoveBugEdges(G,fishy_edges)    
    
    try:
        with pysam.Samfile(bamfile, 'rb') as bam_file:    #once real data, change to 'rb', simulated files are on SAM format
            AddEdges(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len)
    except ValueError:
        with pysam.Samfile(bamfile, 'r') as bam_file:
            AddEdges(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len)


#for edge in G.edges():
    #    if G[edge[0]][edge[1]]['nr_reads']:
    #        print G[edge[0]][edge[1]]['gap_dist']
        
    #print G.edges(data=True)            
    return(G,Contigs,Scaffolds,F,scaffold_indexer)

def RemoveBugEdges(G,fishy_edges):
    edges_removed = 0
    for edge_tuple,nr_links in fishy_edges.items():             
        if edge_tuple[1] in G and edge_tuple[0] in G[edge_tuple[1]]:
            if nr_links >= G[edge_tuple[0]][edge_tuple[1]]['nr_links']:
                G.remove_edge(edge_tuple[0],edge_tuple[1])  
                edges_removed += 1 
#print 'Number of BWA buggy edges removed: ', edges_removed           
    return()

def CheckDir(cont_obj1,cont_obj2,alignedread):
    (read_dir,mate_dir) = (not alignedread.is_reverse,not alignedread.mate_is_reverse )
    cont_dir1 = cont_obj1.direction  #if pos : L if neg: R
    #position2 cont2/scaf2                        
    cont_dir2 = cont_obj2.direction
    (gap,scaf_side1,scaf_side2) = PosDirCalculatorPE(cont_dir1,read_dir,0,0,0,0,cont_dir2,mate_dir,0,0,0,0,0)    
    return(scaf_side1,scaf_side2)

def PosDirCalculatorPE(cont_dir1,read_dir,cont1pos,readpos,s1len,cont1_len,cont_dir2,mate_dir,cont2pos,matepos,s2len,cont2_len,read_len):
    #calculates the positions of the two reads on theis respective contig.
    #o_1^i is denoted gap1 in the code (a bit misleading)
    #o_2^i is denoted gap2 in the code (a bit misleading)
    # o^i is then o_1^i+o_2^i and is denoted gap here
    if cont_dir1 and read_dir:
        gap1=s1len-cont1pos-readpos
        read_side1='R'
    if cont_dir2 and mate_dir:
        gap2=s2len-cont2pos-matepos
        read_side2='R'
    if (not cont_dir1) and read_dir:
        gap1=cont1pos+(cont1_len-readpos)
        read_side1='L'
    if (not cont_dir2) and mate_dir:
        gap2=cont2pos+(cont2_len-matepos)
        read_side2='L'
    if cont_dir1 and not read_dir:
        gap1=cont1pos + readpos + read_len
        read_side1='L'
    if cont_dir2 and not mate_dir:
        gap2=cont2pos + matepos + read_len
        read_side2='L'
    if not cont_dir1 and not read_dir:
        gap1= s1len - cont1pos - (cont1_len-readpos -read_len)
        read_side1='R'
    if not cont_dir2 and not mate_dir:
        gap2= s2len - cont2pos - (cont2_len-matepos -read_len)
        read_side2='R'
    obs=gap1+gap2
    if read_side1 == 'L':
        scaf_side1 = 'L'
    if read_side2 == 'L':
        scaf_side2 = 'L'
    if read_side1 == 'R':
        scaf_side1 = 'R'
    if read_side2 == 'R':
        scaf_side2 = 'R'
    return(obs,scaf_side1,scaf_side2)







