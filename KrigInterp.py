# -*- coding: utf-8 -*-
"""
HBV implementation for BRUE in UK
First oart only considers implementation of precipitation interpolation
"""
import pyximport
pyximport.install()
import VariogramFit

import csv
#import HBV960_0
#import Kriging0_0NZ
import numpy
from numpy import linalg
import cPickle
import random
from pyOpt import ALHSO, Optimization

## get data from csv files
NameSt = []
Loc = []

## location of gauges
with open('SiteInfo.csv','rb') as SiteInfo:
    Lines = csv.reader(SiteInfo)
    Lines.next()
    for row in Lines:
        NameSt.append(row[0])
        Loc.append([float(row[1]),float(row[2])])

print('Gauge Location, Imported')
print('') 

POI = []
CN = []
with open('XY_Basins.csv','rb') as POIf:
    Lines = csv.reader(POIf)
    Lines.next()
    for row in Lines:
        POI.append([float(row[0]),float(row[1])])
        CN.append(int(row[2]))
    
print('Sampling Points, Imported')
print('')    
        
Prec = []
Flow = []
Temper = []
LTTemp = []
Et = []

with open('BrueMauri.csv','rb') as Data:
    Lines = csv.reader(Data)
    Lines.next()
    for row in Lines:
        Prec.append([float(x) for x in row[8:]])
print('Data, Imported')
print('')

## Removal of no precipitation events
WetMeas = []
for i in xrange(0,len(Prec)):
    if numpy.max(Prec[i]) > 3:
        WetMeas.append(Prec[i])    

## Measurement covariance
CovMea = numpy.cov(numpy.transpose(WetMeas))
print('Determinant Covariance Matrix', str(linalg.det(CovMea)))
print('')

#Semivariogram fit

x0 = [1.0,10.0,1.0,1.0,1.0]

    # Boundaries
Sb = (0.01,400) # Limit for the sill
Rb = (2,20) # Limit for the range
Nb = (0,400) # Limit for the Nugget effect
ab = (0,2) # Limit for a in power variogram
vb = (0,1000) # Limit for Matern v parameters

VecBound = (Sb,Rb,Nb,ab,vb)

## Kriging interpolation Precipitation
Z = []
ZAvg = []
SP = []
f = numpy.zeros(len(Prec[0]))

#def PreKrig (Loc,CovMea,Sb,Rb,Nb,ab,vb):
Dis = numpy.zeros((len(Loc),len(Loc)))
for i in xrange(0, len(Loc)):
    j = 0
    for j in xrange(0, len(Loc)):
        Dis[i][j] = numpy.sqrt((Loc[i][0]-Loc[j][0])**2 +
                              (Loc[i][1]-Loc[j][1])**2)
print ""
print "Distance Matrix - Done"

## Do Experimental Semivariogram
## [lag, Cov]
SVExp = []
CovPlt = []
LagPlt = []
for i in xrange(0,len(CovMea)-1):
    for j in xrange(i+1,len(CovMea)):
        Cov = CovMea[i][j]
        Lag = Dis[i][j]
        SVExp.append([Lag,Cov])
        CovPlt.append(Cov)
        LagPlt.append(Lag)
        
print ""
print "Experimental semivariogram - Done"

## Set theoretical variogram vector
## Variogram function array
VarFunArr = [VariogramFit.SVExponential, VariogramFit.SVGaussian, 
             VariogramFit.SVSpherical, VariogramFit.SVCubic,
             VariogramFit.SVPentaspherical, VariogramFit.SVSinehole, 
             VariogramFit.SVPower, VariogramFit.SVMatern]
              
## Variogram Optimisation function array
#    optFunArr = [optSVExponential, optSVGaussian, optSVSpherical, optSVCubic,
#                  optSVPentaspherical, optSVSinehole, optSVPower, optSVMatern]

## names of functions
optFunNam = ['Exponential','Gaussian','Spherical','Cubic',
             'Pentaspherical','Sinehole','Power','Matern']

print ''
print 'Initialising Variogram fit'

sr = random.uniform(Sb[0],Sb[1])
rr = random.uniform(Rb[0],Rb[1])
nr = random.uniform(Nb[0],Nb[1])
ar = random.uniform(ab[0],ab[1])
vr = random.uniform(vb[0],vb[1])

Var = []
Res = []
Mdl = [] 
#    x0 = [sr,rr,nr,ar,vr] #Sill, Range, Nugget, Rank, v exponent
## optimisation begin

def OptFun(x,*args):
#    print len(SVExp)
    F, g, fail = VariogramFit.optFunMaster(x,SVExp,j,VarFunArr)
    if F == 9999:
        fail = 1
    else:
        Var.append(x)
        Res.append(F)
        Mdl.append(j)
    return F, g, fail
    
for j in xrange(0,len(VarFunArr)):    

    VarProb = Optimization('Variogram Fitting: ' + optFunNam[j], OptFun)
    VarProb.addObj('RMSE')
    VarProb.addVar('Sill','c',lower=Sb[0],upper=Sb[1],value=sr)
    VarProb.addVar('Range','c',lower=Rb[0],upper=Rb[1],value=rr)
    VarProb.addVar('Nugget','c',lower=Nb[0],upper=Nb[1],value=nr)
    VarProb.addVar('Exponent a','c',lower=ab[0],upper=ab[1],value=ar)
    VarProb.addVar('Rank v','c',lower=vb[0],upper=vb[1],value=vr)
    
    print ''
    print 'Variogram Fitting ' + optFunNam[j]
    
    args = (SVExp, j, VarFunArr, Var, Res, Mdl)
    
    optmz = ALHSO()
    optmz(VarProb)
    print ''
    print VarProb.solution(0)
    
    ## Best of each semivariograms

## Get position of best semivariogram
k = numpy.argmin(Res)
xopt = Var[k]
ModOpt = Mdl[k]
del Var
del Res
del Mdl

print('')
print("Semivariogram - Done!")
#return VarFunArr,xopt,k,ModOpt


#VarFunArr,xopt,k,ModOpt = PreKrig(Loc,CovMea,Sb,Rb,Nb,ab,vb)
CatNum = []
for i in xrange(max(CN),max(CN)+1):
    POIC = []
    Z = []
    ZAvg = []
    SP = []
#    Z.append('Catchment '+str(i))
#    SP.append('Catchment '+str(i))
#    ZAvg.append('Catchment '+str(i))
    
    for j in xrange(0,len(POI)):
        if i == CN[j]:
            POIC.append(POI[j])

    for ii in xrange(0,len(Prec)):
        if max(Prec[ii]) == 0:
            Z.append(f)
            SP.append(f)
            ZAvg.append(0)
    #        print 'ok, sacando el 0'
        else:
            TempRes = VariogramFit.KrigInterp(ModOpt,POIC,Loc,VarFunArr,xopt,k,CovMea,Prec[ii])
            Z.append(TempRes[0])
            SP.append(TempRes[1])
            ZAvg.append(numpy.average(Z[ii]))
#            CatNum.append(i)
    #        print 'ok, sacando interpolacion'
        if ii%100 == 0:
            print 'Interpolated precipitation register '+ str(ii)
    
    print 'Next Catchment'
    

    with open('AvgPrec-cat'+str(i)+'.pkl','w') as DataFile:
        cPickle.dump(ZAvg,DataFile)
    print 'Average Precipitation Pickled as AvgPrec.pkl'
    
    with open('InterpErr-cat'+str(i)+'.pkl','w') as DataFile:
        cPickle.dump(SP,DataFile)
    print 'Interpolation error pickled as InterpErr.pkl'

