# -*- coding: utf-8 -*-
"""
Created on Thu Jun 06 15:15:57 2013

@author: chaco3
"""

import math
import numpy
from scipy import special
#from scipy.stats import gamma
from numpy import linalg



## Performance Metrics

def RMSE(x,y):
    Erro = numpy.square(numpy.subtract(x,y))
    if Erro.any < 0:
        return 9999
    cdef float F = numpy.sqrt(1.*sum(Erro)/len(x))
    return F

## Semivariogram lists
    # h = Lag dostance
    # S = Partial sill
    # R = Range
    # N = Nugget Effect
    # a = Dimensionless exponent
    # v = Matern parameter
def SVExponential(h,x):
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]   
    cdef float SV
    if h == 0:
        return 0
    if h/R > 3:
        return N+S
    try:
        SV = S - (N + S * (1-numpy.exp(-1.*h/R)))
    except OverflowError:
        SV = N+S
    return SV
    
def SVGaussian(h,x):
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return 0
    if h/R > 1:
        return N+S
    try:
        SV = S - (N + S * (1-numpy.exp(-1.*numpy.square(h)/R)))
    except OverflowError:
        SV = N+S
    return SV

def SVPower(h,x):
    cdef float S = x[0]
    cdef float N = x[2]
    cdef float a = x[3]
    cdef float SV
    if h == 0:
        return 0
    if a > 2:
        a = 2
    try:
        SV = S - (N + S*numpy.power(h,a))
    except OverflowError:
        SV = N+S
    return SV
    
def SVSpherical(h,x):
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV, hr
    if h == 0:
        return 0
    if h > R:
        return N+S
    hr = 1.*h/R
    SV = S - (N + (S * ((3./2)*hr - ((1./2)*(hr**3.)))))
    return SV
 
def SVCubic(h,x):
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return 0
    if h > R:
        return N+S
    SV = S - (N + S * (7*numpy.power((1.*h/R),2.)-(35./4)*numpy.power(1.*h/R,3.)+
        (7./2)*numpy.power(h/R,5)-(3./4)*numpy.power(h/R,7)))
    return SV 

def SVPentaspherical(h,x):
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return 0
    if h > R:
        return N+S
    SV = S - (N + S * ((15./8)*(1.*h/R)-(5./4)*numpy.power(1.*h/R,3)+
        (3./8)*numpy.power(1.*h/R,5)))
    return SV
    
def SVSinehole(h,x):
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return 0
    SV = S - (N + S * (1.-numpy.sin(numpy.pi*1.*h/R)/(numpy.pi*1.*h/R)))
    return SV

def SVMatern(h,x):  
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float v = x[4]
    cdef float SV
    if h == 0:
        return 0
    if v < 0:
        return 9999
    if h/R > 3:
        return N+S
    SV = S - (N + S * (1.-(2./special.gamma(v))*numpy.power((1.*h*numpy.sqrt(v)/R),v)*
        special.kv(2.*h*numpy.sqrt(v)/R,v)))
    return SV
    
def optFunMaster(x,SVExp,j,VarFunArr):
    temp = []
    temp2 = []
    cdef int fail = 0
    for i in xrange (0,len(SVExp)):
        temp.append(VarFunArr[j](SVExp[i][0],x))
        temp2.append(SVExp[i][1])
    cdef float F = RMSE(temp,temp2)
    return F, [], fail

## Serious stuff starts here!

    
def KrigInterp(ModOpt,POI,Loc,VarFunArr,xopt,k,CovMea,MeaRow):
    # Calculate distances between stations        
    SVr = numpy.array(CovMea)
    if linalg.det(CovMea) == 0:
        print('Non-singular covriance matrix - Sorry, cannot invert')
        return(9999,9999)
    InvSVr = linalg.inv(SVr)
    
    POID = numpy.zeros([len(POI),len(Loc)])
    SVm = numpy.zeros([len(POI),len(Loc)])
    # Set interpolation targets
    #POID = [[1487.,0806.,1746.,1044.,1170.,1513.]]
    
    for i in xrange(0,len(POI)):
        for j in xrange(0,len(Loc)):
            # Calculate distance from target to stations
            POID[i][j] = numpy.sqrt((POI[i][0]-Loc[j][0])**2 +
                                   (POI[i][1]-Loc[j][1])**2 )
            # Calculate variance with between observation and targets
            SVm[i][j] = VarFunArr[ModOpt](POID[i][j],xopt)

    WM = []
    Z = []
    SP = []

    for i in xrange(0,len(SVm)):
        WM.append(numpy.dot(InvSVr,numpy.transpose(numpy.array(SVm[i]))))
        Ztemp = numpy.dot(WM[i],MeaRow)
        Z.append(numpy.clip(Ztemp,0,max(MeaRow))) ## cutoff at 0 and max prec
        S = xopt[0]*numpy.ones(len(SVm[i]))-SVm[i]
        SP.append(numpy.dot(WM[i],numpy.transpose(S)))
        

#    
    return Z,SP