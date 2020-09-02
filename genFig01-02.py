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

colors = sns.cubehelix_palette(10) #02
poscol = 6

# General definitions
maxX = 600; startV = 1
xtime = np.linspace(startV, maxX + startV -1,maxX)
xdate = pd.to_datetime(pd.Series(pd.date_range('20200315', periods=maxX)))
rawData = '20200525_Start1503.csv'
randData = '20200526_LAST/'

nSim = 3000

plotNames = ['(a)','(b)','(c)']
#plotNamesPOS = [0.05,0.05,0.05]
#plotIDPOS = [0.5,0.5,0.5]
plt.rc('figure', figsize=(10, 3))
mxCOVID19 = COVID19('SEIRfull', rawData)
mxCOVID19.readStat(randData)
fig, axs = plt.subplots(1, 3,gridspec_kw={'wspace': 0.2})
axins = []
useconf = True; conf = 0.95
betaeta = [[],[],[]]
for indx,state in enumerate(states):
    # Getting the scatter data
    data = mxCOVID19.getFitData(state)
    # Computing the min and max y limits for the time range
    ymin = np.zeros(maxX)+1e8; ymax = np.zeros(maxX)
    minR0 = 1e5; maxR0 = 0
    cnt = 0; remove = 0
    # Analysis of the stat data
    mxCOVID19.getStatData(state)
    statData = mxCOVID19.actStatData['eta'].values
    betaeta[indx].append(mxCOVID19.actStatData['beta'].values)
    betaeta[indx].append(mxCOVID19.actStatData['eta'].values)
    limitsR0 = st.t.interval(conf, len(statData)-1, loc=np.mean(statData), scale=st.sem(statData))
    statData = mxCOVID19.actStatData['beta/eta'].values
    print(conf*100,'% on eta of '+state+':',limitsR0,'. Min R0=',np.min(statData),', Max R0=',np.max(statData),' useconf=', useconf)
    allPars = mxCOVID19.actStatData.copy()
    if (not useconf):
        remove = int(nSim*0.025)
        allPars=allPars.head(-remove) # removing the last n rows
        allPars=allPars.tail(-remove) # removing the first n rows
    for i in range(len(allPars.index)):
        statpars = allPars.iloc[i].values
        #print(statpars[3],limitsR0)
        #if ((not useconf) or (statpars[3]>=limitsR0[0]) and (statpars[3]<=limitsR0[1])):
        if ((not useconf) or (statpars[2]>=limitsR0[0]) and (statpars[2]<=limitsR0[1])):
            minR0 = minR0 if (minR0<statpars[3]) else statpars[3]
            maxR0 = maxR0 if (maxR0>statpars[3]) else statpars[3]
            # Pars -> beta, tau, q, delta, eta, epsilon
            if etaFit:
                # Data input -> id, beta, eta, beta/eta, error, tauOVq, delta, q, epsilon
                mxCOVID19.parVals=tuple([statpars[1],statpars[7]*statpars[5],statpars[7],statpars[6],statpars[2],statpars[8]])
            else:
                # Data input -> id, beta, q, beta/eta, error, tauOVq, delta, eta, epsilon
                mxCOVID19.parVals=tuple([statpars[1],statpars[2]*statpars[5],statpars[2],statpars[6],statpars[7],statpars[8]])
            yVals = mxCOVID19.getModel(xtime)[:,1]
            # Running over x to update max and min values
            cnt+=1
            for xindx in range(maxX):
                ymin[xindx] = ymin[xindx] if (ymin[xindx]<yVals[xindx]) else yVals[xindx]
                ymax[xindx] = ymax[xindx] if (ymax[xindx]>yVals[xindx]) else yVals[xindx]
    print(conf*100,'% on eta of '+state+':',limitsR0,'. Min R0=',minR0,', Max R0=',maxR0,' useconf=', useconf)
    print(state,cnt)
    # Getting the best fitted curve under the given hypothesis
    # Setting pars
    mxCOVID19.parVals=tuple(pars[indx])
    yValsModel = mxCOVID19.getModel(xtime)[:,1]
    # The main axis
    # Plotting the best fit
    axs[indx].plot(xdate,yValsModel,linewidth=0.8,color=colors[poscol])
    # Plotting the shadow probabilities
    axs[indx].fill_between(xdate, ymin, ymax, alpha=0.3,color=colors[poscol])
    axs[indx].tick_params(labelrotation=22)
    axs[indx].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axs[indx].text(0.05,0.95, plotNames[indx], transform=axs[indx].transAxes)
    axs[indx].text(0.45,0.95, state, transform=axs[indx].transAxes)#,color='blue', transform=ax.transAxes)

    # The inset axis
    axins.append(inset_axes(axs[indx], width="30%", height="40%" ,borderpad=1.1))#, loc='lower left', bbox_to_anchor=(0.8, 0.8, 1, 1)))
    # Plotting the data
    axins[-1].scatter(mxCOVID19.realDays.values,data['Infected'],color='gray',s=2)
    # Plotting the shadow probabilities
    axins[-1].fill_between(xtime, ymin, ymax, alpha=0.3,color=colors[poscol])
    # Plotting the best fit
    axins[-1].plot(xtime,yValsModel,linewidth=0.8,color=colors[poscol])
    axins[-1].set_xlim(0,mxCOVID19.realDays.values[-1]) #self.fitData.iloc[0].values
    axins[-1].set_ylim(0,data['Infected'].values[-1]) #self.fitData.iloc[0].values
    axins[-1].tick_params(direction='out', labelsize = 6)
    axins[-1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    print('Maximum for the most optimistic scenario: ', np.max(ymin), " achieved on ", xdate.iloc[np.argmax(ymin, axis=0)])
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.ylabel('Infected')
plt.tight_layout()
plt.savefig('Figure02.png',dpi=600)
plt.show()


smpN=3000
# Details are in the following order: mu,sigma,shift,step,distType
distDetails = [[60,5,0,5,'norm'], # C
               [2,1.4,9,1,'lognorm'], # delta
               [2,1.1,2,1,'lognorm'], # eta
               [0.2,0.03,0,0.05,'norm'], # epsilon
               [10,0.6,0,1,'norm']] #fE0
distPOS = {'C':0,'delta':1,'eta':2,'epsilon':3,'fE0':4, None : 0}
distORDER = [['C','epsilon','fE0'],['eta','delta',None]]
lims = [[[40,80],[0.05,0.35],[7,13]],[[2,10],[8,21]]]

distDetails = [[60,3,0,5,'norm'], # C
               [2,1.4,9,1,'lognorm'], # delta
               [1,0.03,0,0.05,'norm'], # q
               [0.2,0.03,0,0.05,'norm'], # epsilon
               [10,0.6,0,1,'norm']] #fE0
distPOS = {'C':0,'delta':1,'q':2,'epsilon':3,'fE0':4, None : 0}
distORDER = [['C','delta','q'],['fE0','epsilon',None]]
lims = [[[40,80],[8,20],[0.85,1.15]],[[7,13],[0.05,0.35]]]

#r'$\mu=60$\n$\sigma=3$\nType=Normal'
descr = [[r'$C_S$',r'$\delta$',r'$q$'],[r'$f$',r'$\varepsilon$','']]

mxCOVID19 = COVID19('SEIRq', '20200525_Start1503.csv')
dists = mxCOVID19.getRandPars(distDetails,smpN)

plotNames = [['(a)','(b)','(c)'],['(d)','(e)',None]]
plt.rc('figure', figsize=(10, 5))
fig, axs = plt.subplots(2, 3)#,gridspec_kw={'wspace': 0.2, 'hspace' : 0.08},sharex=True)
colors = sns.cubehelix_palette(10)
lines = [0 for i in range(7)]
for irow,row in enumerate(distORDER):
    for icol,distD in enumerate(row):
        if (distD is not None):
            [mu,sigma,shift,step,distType] = distDetails[distPOS[distD]]
            data = dists[distPOS[distD]]
            axs[irow,icol].hist(data, density=True, bins=[v*step+step/2 for v in range(30)], alpha=0.6, color=colors[3], ec=colors[5])
            axs[irow,icol].axvline(mu+shift,color=colors[5])
            axs[irow,icol].set_xlim(lims[irow][icol][0], lims[irow][icol][1])
            #lines[col], = axs[scindex,indx].plot(xdate,yValsModelact,color=colors[col+1]) 
            #lines[-1], = axs[scindex,indx].plot(xdate,yValsModelIni,'-.',color=sns.color_palette()[0])
            #if (indx == 2):
            #    axs[scindex,indx].legend(lines, labels=['Jun','Jul','Aug','Sep','Oct','Nov','Dic','SC'], loc="center right", bbox_to_anchor=(1.4, 0.5),frameon=False)

            #axs[irow,icol].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            axs[irow,icol].text(0.05,0.92, plotNames[irow][icol],transform=axs[irow,icol].transAxes)
            psx = 0.67; deltay = 0.07; sz = 6
            axs[irow,icol].text(psx,0.92-0*deltay, descr[irow][icol],transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(psx,0.92-1*deltay, r'$\mu$='+'{:.2f}'.format(mu+shift),transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(psx,0.92-2*deltay, r'$\sigma$='+'{:.2f}'.format(sigma),transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(psx,0.92-3*deltay, 'Type: '+distType,transform=axs[irow,icol].transAxes,size = sz)
            #axs[irow,icol].tick_params(labelrotation=30,bottom=False)
        else:
            axs[irow,icol].spines['bottom'].set_visible(False) # Hide the bottom axis
            axs[irow,icol].spines['top'].set_visible(False) # Hide the bottom axis
            axs[irow,icol].spines['left'].set_visible(False) # Hide the bottom axis
            axs[irow,icol].spines['right'].set_visible(False) # Hide the bottom axis
            axs[irow,icol].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
            for mv in range(3):
                box = axs[irow,icol-mv].get_position()
                half = (box.x1-box.x0)/2
                box.x0 = box.x0 + half
                box.x1 = box.x1 + half
                axs[irow,icol-mv].set_position(box)

# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

#plt.tight_layout()
plt.savefig('Figure01.png',dpi=600)
plt.show()

for indx in range(3):
    for ith in range(2):
        plt.hist(betaeta[indx][ith], density=True)#, bins=[v*step+step/2 for v in range(30)], alpha=0.6, color=colors[3], ec=colors[5])
        plt.show()