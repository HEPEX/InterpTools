# -*- coding: utf-8 -*-
"""
From pickle to CSV
"""

import csv
import cPickle
#import numpy.core.multiarray as multiarray

for i in xrange(1,69):
    with open('AvgPrec-cat'+str(i)+'.pkl','r') as DataFile:
        Data = cPickle.load(DataFile)
        
    with open('AvgPrec2time-cat'+str(i)+'.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar=',', 
                                quoting=csv.QUOTE_MINIMAL)        
        for j in xrange(len(Data)):
            spamwriter.writerow([Data[j]])



    