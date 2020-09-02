from importCOVID19 import COVID19
import numpy as np

states = [['CMX',  9],['MEX', 15],['JAL', 14],['PUE', 21],['NLE', 19],
          ['COA',  5],['YUC', 31],['TAB', 27],['ROO', 23],['GUA', 11],
          ['AGU',  1],['BCN',  2],['SIN', 25],['SLP', 24],['QUE', 22],
          ['VER', 30],['MIC', 16],['OAX', 20],['HID', 13],['SON', 26],
          ['BCS',  3],['GRO', 12],['CHP',  7],['CHH', 8],['TAM', 28],['MOR', 17], 
          ['NAY', 18],['DUR', 10],['ZAC', 32],['CAM',  4],['TLA', 29],['COL',  6 ]]

statesLST = ['TOT','GR0','GR1','GR2']
statesLST.extend([ st[0] for st in states ])
mxCOVID19 = COVID19('SEIReta', '20200525_Start1503.csv')
mxCOVID19.doStat(20,statesLST,2,'rand',
                 [[60,3,0,5,'norm'],[2,1.4,9,1,'lognorm'],[1,0.03,0,0.05,'norm'],
                  [0.2,0.03,0,0.05,'norm'],[10,0.6,0,1,'norm']],1)
