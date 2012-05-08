import sys
import Contig,Scaffold
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import pi
def GapEstimator(G,Contigs,Scaffolds,mean,sigma,read_length,dValuesTable,edge_support):
    gap_mean=0
    gap_obs=0
    gap_obs2=0
    gap_counter=0
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            c1=edge[0][0]
            c2=edge[1][0]
            c1_len=G.node[edge[0]]['length']
            c2_len=G.node[edge[1]]['length']
            obs_list=G.edge[edge[0]][edge[1]]['gap_dist']
            nr_links=G.edge[edge[0]][edge[1]]['nr_links']
            if nr_links >= edge_support:
                # if long contigs, use pre calculated d_ML fcn
                if (c1_len > mean+6*sigma) and (c1_len > mean+6*sigma):
                    d_ML,stdErr=CalcGapLongContigs(obs_list,nr_links,mean,sigma,read_length,c1_len,c2_len,dValuesTable)
                    #print 'Longa contigs'

                #calculate d_ML for particular contig lengths
                else:
                    d_ML,stdErr=CalcMLvaluesOfdGeneral(obs_list,mean,sigma,read_length,c1_len,c2_len,nr_links)               
                
                gap_obs+=d_ML
                gap_obs2+=d_ML**2
                gap_counter+=1
                print "gap distance: ",d_ML,'between: ', Scaffolds[c1].contigs[0].name,' and ', Scaffolds[c2].contigs[0].name, 'Number of links: ',nr_links #, 'Data obs: ', (nr_links*mean -int(sum(obs_list)))/float(nr_links)
                #,'std_err: ',stdErr,c1_len,c2_len, 'OBS= ',(nr_links*mean -int(sum(obs_list)))/float(nr_links), nr_links
    gap_mean=gap_obs/float(gap_counter)
    gap_sigma=gap_obs2-gap_counter*gap_mean**2
    print "Avg gap: ",gap_mean ,"Std_gap_estimate: ", (gap_sigma/(gap_counter-1))**(0.5) , "nr of gaps: ",gap_counter




def CalcGapLongContigs(obs_list,nr_links,mean,stdDev,readLen,c1_len,c2_len,dValuesTable):  
    #get observation    
    data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    MLgap=dValuesTable[int(round(data_observation,0))][0]
    stdErr=(1.0/(nr_links*dValuesTable[int(round(data_observation,0))][1]))**0.5
    return(MLgap,stdErr)

def PreCalcMLvaluesOfdLongContigs(mean,stdDev,readLen):
    def Nom(z,mean,stdDev):
        nom=-(1+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        return nom
    def Denom(d,readLen,mean,stdDev):
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*(1+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))
        second=stdDev*2.718**(-((mean-d-2*readLen+1)**2)/(float(2*stdDev**2)))
        denom=first+second
        return denom

    def CalcAofd(d):
        #transform to z ~N(0,1)
        z=((d+2*readLen)-mean)/float(stdDev)    
        nom=Nom(z,mean,stdDev)
        denom=Denom(d,readLen,mean,stdDev)
        val=nom/denom              
        return val
    #we cannot 
    d_upper=mean+6*stdDev
    d_lower=-2*stdDev
    dValuesTable={}
    for d in range(d_lower,d_upper+1):
        Aofd=CalcAofd(d)
        obs=d+Aofd*stdDev**2
        obs=int(round(obs,0))
        #print d,obs
        Iofd=1.0/(stdDev**2)+Aofd**2
        #print Iofd,stdDev
        # a table with observation as key and the information function I(d) for d
        # to get std err and CI of ML estimate for d, do calc: sqrt(1.0/(n*I(d)))
        dValuesTable[obs]=(d,Iofd)
        
    return dValuesTable

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
        second=stdDev*( 2.718**(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) - 2.718**(-( (c_min+d+readLen-mean)**2)/(float(2*stdDev**2)))) 
        term1=first+second

        first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev)))  )
        #print 'First: ',first
        second=stdDev*( 2.718**(-( (c_min+c_max+d-mean)**2)/(float(2*stdDev**2))) - 2.718**(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))))
        #print 'Second: ',second
        term3=first+second
        denom=term1+term2+term3
        #print term1,term2,term3
        return denom
    
    denominator=Denominator(d,c1Len,c2Len,readLen)
    nominator,case_nomin=Nominatortest(d,c_min,c_max,c1Len,c2Len,readLen)
    ## Aofdlist=[]
    ## for obs in obs_list:
    ##     nominator,case_nomin=Nominator(obs,c_min,c_max,c1Len,c2Len,readLen)
    ##     #sys.stdout.write(str(case_nomin))
    ##     Aofdlist.append(nominator/float(denominator))
    ## #print 'Nom: ',nominator
    ## Aofd=float(sum(Aofdlist))/len(Aofdlist)
    Aofd=nominator/denominator
    func_of_d=d+Aofd*stdDev**2
    return func_of_d

def CalcMLvaluesOfdGeneral(obs_list,mean,stdDev,readLen,c1Len,c2Len,nr_links):
    #get observation    
    data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    #do binary search among values
    d_upper=mean+6*stdDev
    d_lower=-2*stdDev
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
