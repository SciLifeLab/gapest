import sys
import Contig,Scaffold
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import pi
from math import exp
def GapEstimator(G,Contigs,Scaffolds,mean,sigma,read_length,edge_support):
    gap_mean=0
    gap_obs=0
    gap_obs2=0
    gap_counter=0
    print 'CONTIG1\tCONTIG2\tGAP_ESTIMATION\tNUMBER_OF_OBSERVATIONS\tWARNINGS/ERRORS'
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            c1=edge[0][0]
            c2=edge[1][0]
            c1_len=G.node[edge[0]]['length']
            c2_len=G.node[edge[1]]['length']
            obs_list=G.edge[edge[0]][edge[1]]['gap_dist']
            nr_links=G.edge[edge[0]][edge[1]]['nr_links']
            #pre check for large deviations in obs list
            sorted_observations = sorted(obs_list)
            smallest_obs_mean = sum(sorted_observations[0:10])/10.0
            largest_obs_mean = sum(sorted_observations[-10:])/10.0
            #print largest_obs_mean,smallest_obs_mean
            if nr_links >= edge_support and largest_obs_mean-smallest_obs_mean < 6*sigma:
                d_ML,stdErr=CalcMLvaluesOfdGeneral(obs_list,mean,sigma,read_length,c1_len,c2_len,nr_links)               
                
                gap_obs+=d_ML
                gap_obs2+=d_ML**2
                gap_counter+=1
                warn = 0
                if c1_len < 2*sigma and c2_len < 2*sigma:
                    warn = 1
                if warn == 1:
                    print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t'+str(d_ML) +'\t'+str(nr_links)+'\tw'
                else:
                    print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t'+str(d_ML) +'\t'+str(nr_links)+'\t-'


            elif nr_links < edge_support:
                print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t.\t.\te1'
            elif largest_obs_mean-smallest_obs_mean < 6*sigma:
                print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t.\t.\te2'
                    
    print 'w : Both contig lengths were smaller than 2*std_dev of lib (heuristic threshold set by me from experience). This can give shaky estimations in some cases.'
    print 'e1 : No gap was calculated. Number of links were lower than specified min nr of links parameter: -e <min nr links> (default 10). Lower this value if necessary (estimations may become unstable)'
    print 'e2 : No gap was calculated. The spread of the links throughout the contig is to far (i.e largest_obs_mean-smallest_obs_mean < 6*sigma ), suggesting errors in mapping on this region.',
                
    gap_mean=gap_obs/float(gap_counter)
    gap_sigma=gap_obs2-gap_counter*gap_mean**2

def funcDGeneral(obs_list,d,mean,stdDev,c1Len,c2Len,readLen):
    #get observation    
    #data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    c_min=min(c1Len,c2Len)
    c_max=max(c1Len,c2Len)

    def Nominatortest(d,c_min,c_max,c1Len,c2Len,readLen):
        nomin1=( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) +erf((mean-d-c_max-readLen-1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        nomin2=-(erf((c_min+d+readLen+1-mean)/(2**0.5*float(stdDev)))+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        nomin=nomin1+nomin2
        case_nomin=0
        return nomin,case_nomin
    
    def Denominator(d,c1Len,c2Len,readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term2=(c_min-readLen+1)/2.0*(erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))- erf((c_min+d+readLen-mean)/((2**0.5)*stdDev))   )
    
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*( erf((c_min+d+readLen-mean)/(2**0.5*float(stdDev))) - erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev)))  )
        second=stdDev*( exp(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_min+d+readLen-mean)**2)/(float(2*stdDev**2)))) 
        term1=first+second

        first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev)))  )
        #print 'First: ',first
        second=stdDev*( exp(-( (c_min+c_max+d-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))))
        #print 'Second: ',second
        term3=first+second
        denom=term1+term2+term3
        #print term1,term2,term3
        return denom
    
    denominator=Denominator(d,c1Len,c2Len,readLen)
    nominator,case_nomin=Nominatortest(d,c_min,c_max,c1Len,c2Len,readLen)
    Aofd=nominator/denominator
    func_of_d=d+Aofd*stdDev**2
    return func_of_d

def CalcMLvaluesOfdGeneral(obs_list,mean,stdDev,readLen,c1Len,c2Len,nr_links):
    #get observation    
    data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    #do binary search among values
    d_upper=int(mean+2*stdDev-2*readLen)
    d_lower=-10*stdDev
    while d_upper-d_lower>1:
        d_ML=(d_upper+d_lower)/2.0
        func_of_d=funcDGeneral(obs_list,d_ML,mean,stdDev,c1Len,c2Len,readLen)
        #print d_ML,func_of_d
        if func_of_d>data_observation:
            d_upper=d_ML
        else:
            d_lower=d_ML

    
    d_ML=(d_upper+d_lower)/2.0
    return int(round(d_ML,0)),0
