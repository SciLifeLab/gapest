'''
Created on Sep 23, 2011

@author: ksahlin
'''

import pysam
import Contig,Scaffold
import networkx as nx

def PE(Contigs,Scaffolds,bamfile,mean,scaffold_indexer,F,read_len):
    G=nx.Graph()
    print 'Parsing BAM file...'
    #read_len=50
    #informative_pair={81:(False,True),97:(True,False),113:(False,False),65:(True,True)}
    #I switched to look at mates instead since BWA can give false flag combinations for
    # read-mate when read is mapped but not mate eg 97-149 81-165. But the reverse
    #does not happen.
    informative_pair={161:(True,False),145:(False,True),129:(True,True),177:(False,False)}
    #threshold=800
    with pysam.Samfile(bamfile, 'r') as bam_file:    #once real data, change to 'rb', simulated files are on SAM format
        #Clean contig_library
        singeled_out=0
        cont_lengths= bam_file.lengths
        cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
        #print cont_lengths
        cont_names = bam_file.references
####### WHEN ADDING SHORTER CONTIGS NOT INCLUDED IN THE SCAFFOLDING, 
####### WE NEED TO ALSO INITIALIZE OBJECTS FOR THESE, THIS SHOULD BE DONE SOMEWHERE HERE
        for i in range(0,len(cont_names)):
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
        print 'Nr of contigs/scaffolds included in scaffolding: '+ str(len(Scaffolds))#,Scaffolds.keys()
        
        for scaffold_ in Scaffolds:
            G.add_edge((scaffold_,'L'),(scaffold_,'R'),nr_links=None)    #this is a scaffold object but can be both a single contig or a scaffold.
                                                       
                
        # Create the link edges in the graph by fetching info from bam file

            
        for alignedread in bam_file:
            flag_type=alignedread.flag
            if flag_type in informative_pair:
                contig1=bam_file.getrname(alignedread.rname)
                contig2=bam_file.getrname(alignedread.mrnm)
                if contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
                    (read_dir,mate_dir)=informative_pair[flag_type]
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
                    (gap,scaf_side1,scaf_side2)=PosDirCalculatorPE(cont_dir1,read_dir,cont1_pos,readpos,s1len,cont1_len,cont_dir2,mate_dir,cont2_pos,matepos,s2len,cont2_len,read_len)   
                    if (scaf2,scaf_side2) not in G[(scaf1,scaf_side1)]:
                        G.add_edge((scaf2,scaf_side2),(scaf1,scaf_side1),nr_links=1,gap_dist=[gap])
                        #print 'Added edge'
                    else:
                        G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['nr_links'] += 1
                        #print 'edge'
                        G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['gap_dist'].append(gap)                         
                                
                elif contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
########################Use to validate scaffold herein previous step here
                    pass
    #for edge in G.edges():
    #    if G[edge[0]][edge[1]]['nr_reads']:
    #        print G[edge[0]][edge[1]]['gap_dist']
        
    #print G.edges(data=True)            
    return(G,Contigs,Scaffolds,F,scaffold_indexer)





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
    gap=gap1+gap2
    if read_side1 == 'L':
        scaf_side1 = 'L'
    if read_side2 == 'L':
        scaf_side2 = 'L'
    if read_side1 == 'R':
        scaf_side1 = 'R'
    if read_side2 == 'R':
        scaf_side2 = 'R'
    return(gap,scaf_side1,scaf_side2)







