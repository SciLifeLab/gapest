#import numpy
import sys
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import pi
from numpy import *
import matplotlib.pyplot as plt
#print erf(1),norm.cdf(2)
d=-400
readLen=50
mean=650
stdDev=65
def funcD(d,mean,stdDev):
    #transform to z ~N(0,1)
    z=((d+2*readLen)-mean)/float(stdDev)
    #print 'Z: ',z
    def Nom(z,mean,stdDev):
        #nom=(1-norm.cdf(z))*(2*pi)**0.5 #+ norm.pdf(z)*(2*pi)**0.5
        nom=-(1+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        #print 'Nom: ',nom
        return nom
    def Denom(d,readLen,mean,stdDev):
        #denom=stdDev*norm.pdf(((d+2*readLen)-mean)/float(stdDev)) + (mean+1-d-2*readLen)*(1-norm.cdf((d+2*readLen-mean)/float(stdDev)))
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*(1+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))
        second=stdDev*2.718**(-((mean-d-2*readLen+1)**2)/(float(2*stdDev**2)))
        denom=first+second
        return denom
    
    nom=Nom(z,mean,stdDev)
    denom=Denom(d,readLen,mean,stdDev)
    val=nom/denom              
    return val

def funcDGeneral(d,mean,stdDev,c1Len,c2Len):
    c_min=min(c1Len,c2Len)
    c_max=max(c1Len,c2Len)
    
    #transform to z ~N(0,1)
    #z=((d+2*readLen)-mean)/float(stdDev)

    def Nominatortemp(d,c_min,c_max,c1Len,c2Len,readLen):
        #case 2-3: (c_max + d + r >= x)  ^ (c_min + d + r <= x)
        #if c_max+readLen >= obs and c_min + readLen <=obs:
        #    nomin=0
        #    case_nomin=1
        #elif c_max + readLen < obs and c1Len+c2Len + 2*readLen <= 2*obs:
        nomin1=( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) +erf((mean-c_max-readLen-1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        #case_nomin=2
        #often case 1: paired read restricted by insert size
        #else:
        nomin2=-(erf((c_min+readLen+1-mean)/(2**0.5*float(stdDev)))+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        #    case_nomin=3
        nomin=nomin1+nomin2
        case_nomin=0
        return nomin,case_nomin

    def Nominator(obs,c_min,c_max,c1Len,c2Len,readLen):
        #case 2-3: (c_max + d + r >= x)  ^ (c_min + d + r <= x)
        if c_max+readLen >= obs and c_min + readLen <=obs:
            nomin=0
            case_nomin=1
        elif c_max + readLen < obs and c1Len+c2Len + 2*readLen <= 2*obs:
            nomin=( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) +erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
            case_nomin=2
        #often case 1: paired read restricted by insert size
        else:
            nomin=-(erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev)))+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
            case_nomin=3
        return nomin,case_nomin 
    
    def Denominator(d,c1Len,c2Len,readLen):
        #possibleCases=[obs-2*readLen+1,c_min-readLen+1,c1Len+c2Len-obs+1]
        #case=possibleCases.index(min(possibleCases))
        #denominator=min(possibleCases)
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term2=(c_min-readLen+1)/2.0*(erf((c_max+d+readLen-mean)/(sqrt(2)*stdDev))- erf((c_min+d+readLen-mean)/(sqrt(2)*stdDev))   )
    
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*( erf((c_min+d+readLen-mean)/(2**0.5*float(stdDev))) - erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev)))  )
        second=stdDev*( exp(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_min+d+readLen-mean)**2)/(float(2*stdDev**2)))) 
        term1=first+second

        first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev)))  )
        #print 'First: ',first
        second=stdDev*( exp(-( (c_min+c_max+d+1-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))))
        #print 'Second: ',second
        term3=first+second
        denom=term1+term2+term3
        #print 'First: ', term1,term2,'Third: ',term3, 'x: ',d
        return denom
    
    denominator=Denominator(d,c1Len,c2Len,readLen)
    nominator,case_nomin=Nominatortemp(d,c_min,c_max,c1Len,c2Len,readLen)
    ## Aofdlist=[]
    ## for obs in range(300,400):
    ##     #nominator,case_nomin=Nominator(20,c_min,c_max,c1Len,c2Len,readLen)
    ##     #sys.stdout.write(str(case_nomin))
    ##     Aofdlist.append(nominator/float(denominator))
    ## #print d+2*readLen,'Nom: ',nominator, 'Denom: ',denominator, 'Term', (nominator/denominator)**2
    ## Aofd=float(sum(Aofdlist))/len(Aofdlist)
    Aofd=nominator/denominator
    return Aofd



Aofd=0
y=[]
#x=[]
c1Len=300
c2Len=300
for i in range(0,140):
    #d from 0 to 900
    a_prev=Aofd
    #Aofd1=funcD(d,mean,stdDev)
    Aofd=funcDGeneral(d,mean,stdDev,c1Len,c2Len)
    #Aofd=Aofd1-Aofd2
    #print 'Aofd: ',Aofd*stdDev,d
    y.append(d+Aofd*stdDev**2)
    #x.append(d)
    #if (d+Aofd*stdDev**2 - mean-2*readLen) < 10 and (d+Aofd*stdDev**2 - mean-2*readLen) > -10 :
    #    upper_quantile=(d-mean)/float(stdDev)
    #    print upper_quantile
    print 'func: ', d+Aofd*stdDev**2, 'd: ',d, Aofd #'I(d): n*',((float(1)/stdDev**2) + Aofd**2)  #, (Aofd-a_prev)*stdDev**2
    d=d+10
    
x = arange(-400,1000,10)
plt.plot(x,y,'k-',label='$\mu=650,\ \sigma=300$')
y=[]
for i in range(0,140):
    y.append(1000)
plt.plot(x,y,'r',label='max obs')
y=[]
for i in range(0,140):
    y.append(1100)
plt.plot(x,y,'r',label='max obs')
plt.xlabel('d_ML')
plt.ylabel('observation')
#plt.text(100, .025, r'$\mu=650,\ \sigma=300$')
plt.title('ML function of d')
#plt.plot(x,y,'r',label='max obs')
#plt.legend(loc=0)

