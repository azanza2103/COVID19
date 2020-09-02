from importCOVID19 import COVID19
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rc('figure', figsize=(10, 5))
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# For statistical plots
import seaborn as sns
import scipy.stats as st
from scipy.stats import norm

# Reading parameters fit up to 2020/05/15
fitpars = pd.read_csv('20200525fits.csv')
# Parameters up to 2020/05/12
# First row is for eta=1/4 and second one for eta=1/5 [beta,q]
# Third row is for eta fiting with q=1 [beta,eta] up to 2020/05/15
res=[[[0.744,0.245],[0.734,0.105],[0.744,0.144]],
     [[0.617,0.162],[0.590,0.065],[0.610,0.097]],
     [fitpars[fitpars['ID']=='TOT'].values[0][1:3],
      fitpars[fitpars['ID']=='CMX'].values[0][1:3],
      fitpars[fitpars['ID']=='MEX'].values[0][1:3]]]
etaAct = 1/4; etaFit = True
posInres = 2 if etaFit else (0 if (etaAct==1/4) else 1)
if etaFit:
    q = 1
    # TOT
    # beta, tau, q, delta, eta, epsilon
    parsTOT = [res[posInres][0][0],q*2/3,q,1/10,res[posInres][0][1],0.2]
    # CMX
    # beta, tau, q, delta, eta, epsilon
    parsCMX = [res[posInres][1][0],q*2/3,q,1/10,res[posInres][1][1],0.2]
    # MEX
    # beta, tau, q, delta, eta, epsilon
    parsMEX = [res[posInres][2][0],q*2/3,q,1/10,res[posInres][2][1],0.2]
else:
    # TOT
    # beta, tau, q, delta, eta, epsilon
    parsTOT = [res[posInres][0][0],res[posInres][0][1]*2/3,res[posInres][0][1],1/10,etaAct,0.2]
    # CMX
    # beta, tau, q, delta, eta, epsilon
    parsCMX = [res[posInres][1][0],res[posInres][1][1]*2/3,res[posInres][1][1],1/10,etaAct,0.2]
    # MEX
    # beta, tau, q, delta, eta, epsilon
    parsMEX = [res[posInres][2][0],res[posInres][2][1]*2/3,res[posInres][2][1],1/10,etaAct,0.2]
states = ['TOT','CMX','MEX']
pars = [parsTOT,parsCMX,parsMEX]

# General definitions
maxX = 600; startV = 1
xtime = np.linspace(startV, maxX + startV -1,maxX)
xdate = pd.to_datetime(pd.Series(pd.date_range('20200315', periods=maxX)))
rawData = '20200525_Start1503.csv'
randData = '20200526_LAST/'

nSim = 3000
mxCOVID19rand = COVID19('SEIRfull', rawData)
mxCOVID19rand.readStat(randData)

plotNames = [['(a)','(b)','(c)'],['(d)','(e)','(f)'],['(g)','(h)','(i)']]
plotNamesPOS = [[0.95e7,7e5,1.38e6],[4.5e6,3.6e5,7.5e5],[1.14e6,0.8e5,1.6e5]]
plotIDPOS = [0.95e7,7e5,1.38e6]
plt.rc('figure', figsize=(10, 3))
states = ['TOT','CMX','MEX']
schemes = [6,5,5]
bestStrategy = [5,4,5]
startMonths = [232,171,201]
startMonths = [220,176,210]
startMonths = [240,188,233]
pars = [parsTOT,parsCMX,parsMEX]
fig, axs = plt.subplots(1, 3,gridspec_kw={'wspace': 0.2, 'hspace' : 0.08},sharex=True)

colors = sns.light_palette("navy", 10) #04
colors = sns.color_palette("GnBu_d",10) #01
colors = sns.light_palette("navy", 10,reverse=True) #03
colors = sns.cubehelix_palette(30) #02

twinaxs = [0,0,0]
lines = [0 for i in range(7)]
palete = sns.color_palette("GnBu_d")
useconf = True; conf = 0.95
for indx,state in enumerate(states):
    # The first row plots are for quarantine lifting on June 1rst
    mxCOVID19 = COVID19('SEIRfull', rawData)
    # beta, tau, q, delta, eta, epsilon
    mxCOVID19.parVals=tuple(pars[indx])
    data = mxCOVID19.getFitData(state)
    mxCOVID19.eta = etaAct
    yValsModelIni = mxCOVID19.getModel(xtime)[:,1]
    mxCOVID19.tauInt = True; mxCOVID19.tauLinear = False
    mxCOVID19.tauStart = 79 # 1ro junio
    mxCOVID19.tauStep = 20; mxCOVID19.tauDeltad = 1
    mxCOVID19.eta = etaAct
    yValsModel = mxCOVID19.getModel(xtime)[:,1]
    axs[indx].plot(xdate,yValsModel,color=colors[4]) #9
    axs[indx].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axs[indx].text(0.5,0.92,state,transform=axs[indx].transAxes)
    axs[indx].text(0.05,0.92, plotNames[0][indx],transform=axs[indx].transAxes)
    axs[indx].tick_params(labelrotation=22,bottom=False)

# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.ylabel('Infected')

plt.tight_layout()
plt.savefig('Figure06.png',dpi=600)
plt.show()


plt.rc('figure', figsize=(10, 6))
fig, axs = plt.subplots(2, 3,gridspec_kw={'wspace': 0.2, 'hspace' : 0.08},sharex=True)

colors = sns.light_palette("navy", 10) #04
colors = sns.color_palette("GnBu_d",10) #01
#colors = sns.light_palette("navy", 10,reverse=True) #03
#colors = sns.cubehelix_palette(30) #02
colors = sns.color_palette("Paired",8) #01

twinaxs = [0,0,0]
lines = [0 for i in range(7)]
palete = sns.color_palette("GnBu_d")
useconf = True; conf = 0.95
for indx,state in enumerate(states):
    # The first row plots are for quarantine lifting on June 1rst
    mxCOVID19 = COVID19('SEIRfull', rawData)
    # beta, tau, q, delta, eta, epsilon
    mxCOVID19.parVals=tuple(pars[indx])
    data = mxCOVID19.getFitData(state)
    mxCOVID19.eta = etaAct
    yValsModelIni = mxCOVID19.getModel(xtime)[:,1]
    # The second row are for the best strategy (6 for total and 5 for the other two) with different start months
    #for col,start in enumerate([79,109,140,171,201,232,262]):
    for col,start in enumerate([79,109,140,171,201,232]):
        mxCOVID19 = COVID19('SEIRfullV01', rawData,[schemes[indx],20,start,-1,1])
        # beta, tau, q, delta, eta, epsilon
        mxCOVID19.parVals=tuple(pars[indx])
        data = mxCOVID19.getFitData(state)
        #mxCOVID19.eta = etaAct
        yValsModelact = mxCOVID19.getModel(xtime)[:,1]
        lines[col], = axs[1-1,indx].plot(xdate,yValsModelact,color=colors[col]) 
    lines[-1], = axs[1-1,indx].plot(xdate,yValsModelIni,'-.',color=sns.color_palette()[7])
    if (indx == 2):
        #axs[1,indx].legend(lines, labels=['Jun','Jul','Aug','Sep','Oct','Nov','Dic','SC'], loc="center right", bbox_to_anchor=(1.4, 0.5),frameon=False)
        axs[1-1,indx].legend(lines, labels=['Jun','Jul','Aug','Sep','Oct','Nov','SC'], loc="center right", bbox_to_anchor=(1.4, 0.5),frameon=False)

    axs[1-1,indx].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axs[1-1,indx].text(0.05,0.92, plotNames[1-1][indx],transform=axs[1-1,indx].transAxes)
    axs[1-1,indx].tick_params(labelrotation=30,bottom=False)

    # The third row shows the percentage strategy and the best result
    mxCOVID19 = COVID19('SEIRfullV01', rawData,[schemes[indx],20,startMonths[indx],-1,1])
    mxCOVID19.readStat(randData)
    mxCOVID19.getStatData(state)
    # beta, tau, q, delta, eta, epsilon
    mxCOVID19.parVals=tuple(pars[indx])
    data = mxCOVID19.getFitData(state)
    yValsModelact = mxCOVID19.getModel(xtime)[:,1]
    axs[2-1,indx].plot(xdate,yValsModelact,color=colors[bestStrategy[indx]])
    # Generating the 'shadow' zone
    # Computing the min and max y limits for the time range
    ymin = np.zeros(maxX)+1e8; ymax = np.zeros(maxX)
    minR0 = 1e5; maxR0 = 0
    cnt = 0; remove = 0
    # Analysis of the stat data
    statData = mxCOVID19.actStatData['eta'].values
    limitsR0 = st.t.interval(conf, len(statData)-1, loc=np.mean(statData), scale=st.sem(statData))
    statData = mxCOVID19.actStatData['beta/eta'].values
    print(conf*100,'% on eta of '+state+':',limitsR0,'. Min R0=',np.min(statData),', Max R0=',np.max(statData),' useconf=', useconf)
    allPars = mxCOVID19.actStatData.copy()
    #if (not useconf):
    #    remove = int(nSim*0.025)
    #    allPars=allPars.head(-remove) # removing the last n rows
    #    allPars=allPars.tail(-remove) # removing the first n rows
    for i in range(len(allPars.index)):
        statpars = allPars.iloc[i].values
        #print(statpars[3],limitsR0)
        if ((not useconf) or (statpars[2]>=limitsR0[0]) and (statpars[2]<=limitsR0[1])):
            #print(statpars)
            minR0 = minR0 if (minR0<statpars[3]) else statpars[3]
            maxR0 = maxR0 if (maxR0>statpars[3]) else statpars[3]
            # Pars -> beta, tau, q, delta, eta, epsilon
            mxCOVID19.fE0 = statpars[9]
            if etaFit:
                # Data input -> id, beta, eta, beta/eta, error, tauOVq, delta, q, epsilon
                mxCOVID19.parVals=tuple([statpars[1],statpars[7]*statpars[5],statpars[7],statpars[6],statpars[2],statpars[8]])
            else:
                # Data input -> id, beta, q, beta/eta, error, tauOVq, delta, eta, epsilon
                print('whooooooooo')
                mxCOVID19.parVals=tuple([statpars[1],statpars[2]*statpars[5],statpars[2],statpars[6],statpars[7],statpars[8]])
                
            yVals = mxCOVID19.getModel(xtime)[:,1]
            # Running over x to update max and min values
            cnt+=1
            for xindx in range(maxX):
                ymin[xindx] = ymin[xindx] if (ymin[xindx]<yVals[xindx]) else yVals[xindx]
                ymax[xindx] = ymax[xindx] if (ymax[xindx]>yVals[xindx]) else yVals[xindx]
    print(conf*100,'% on eta of '+state+':',limitsR0,'. Min R0=',minR0,', Max R0=',maxR0,' useconf=', useconf)
    print(state,cnt)
    # Plotting the shadow probabilities
    axs[2-1,indx].fill_between(xdate, ymin, ymax, alpha=0.3,color=colors[bestStrategy[indx]])
    print('Maximum for the most optimistic scenario: ', np.max(ymin), " achieved on ", xdate.iloc[np.argmax(ymin, axis=0)])
    
    
    axs[2-1,indx].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axs[2-1,indx].text(0.05,0.88, plotNames[2-1][indx],transform=axs[2-1,indx].transAxes)
    axs[2-1,indx].tick_params(labelrotation=30)
    twinaxs[indx] = axs[2-1,indx].twinx()  # instantiate a second axes that shares the same x-axis
    mxCOVID19.breakDay = startMonths[indx]
    perVals = [ mxCOVID19.getqtau(val,schemes[indx],-1)[2] for val in xtime ]
    twinaxs[indx].plot(xdate, perVals, color=colors[7],linewidth=0.8, alpha=0.6)
    if (indx<2):
        twinaxs[indx].tick_params(labelcolor='none',right=False)#labelrotation=30)
    else:
        twinaxs[indx].set_ylabel('Social confinement (\%)')
        xmin, xmax, yminv, ymaxv = twinaxs[indx].axis()
        twinaxs[0].set_ylim(yminv,ymaxv)

    #twinaxs[indx].tick_params(axis='y', labelcolor=colors[7])
    #break

# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.ylabel('Infected')

#plt.tight_layout()
plt.savefig('Figure07.png',dpi=600)
plt.show()


