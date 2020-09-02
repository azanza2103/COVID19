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
# Excluded states
toexclude = ['TOT','CMX','MEX','GR0','GR1','GR2','COA','BCS','DUR','ZAC','COL','ROO','AGU']
# Excluding further states for high R0 values
#toexclude.extend(['OAX','SON','GRO','HID','SLP','QUE','CAM']) # From data 15/05
toexclude.extend(['OAX','GRO','SLP','CHH','MOR','BCN','CHP']) # From data 25/05
mainpars = fitpars.loc[(fitpars['ID'] == 'TOT') | (fitpars['ID'] == 'CMX') | (fitpars['ID'] == 'MEX')]
print(mainpars)
fitpars = fitpars[~fitpars['ID'].isin(toexclude)]
# Sorting pars in terms of the maximum infected value
fitpars.sort_values(by=['maxI'], inplace=True,ascending=False)
print('Number of states to plot:', len(fitpars.index))

# Get the parVals for a given fitpars row
def getPars(row,stat=False):
    if stat:
        return (row[1],row[7]*row[5],row[7],row[6],row[2],row[8])
    # row = [ID, beta, eta]
    # beta, tau, q, delta, eta, epsilon
    q = 1; epsilon = 0.2; delta = 1/10
    return (row[1],q*2/3,q,delta,row[2],epsilon)

# General definitions
maxX = 600; startV = 1
xtime = np.linspace(startV, maxX + startV -1,maxX)
xdate = pd.to_datetime(pd.Series(pd.date_range('20200315', periods=maxX)))
rawData = '20200525_Start1503.csv'
randData = '20200526_LAST/'
colors = sns.cubehelix_palette(10) #02
poscol = 6
mxCOVID19 = COVID19('SEIRfull', rawData)
mxCOVID19.readStat(randData)
plt.rc('figure', figsize=(12, 12))
fig, axs = plt.subplots(4, 4, gridspec_kw={'wspace': 0.2, 'hspace' : 0.08}, sharex=True, sharey=True)
useconf = True; conf = 0.95

nstate = 0; axins = []; source = mainpars; maxstates = len(source.index)
outres = []
for irow in range(4):
    for icol in range(4):
        print(irow,icol,nstate,maxstates)
        for ps in range(maxstates+1):
            if (nstate<(len(source.index)+maxstates)):
                if ps==maxstates:
                    source = fitpars
                    maxstates = 0
                print(ps,maxstates)
                row = source.iloc[nstate if maxstates==0 else ps].values
                state = row[0]
                print('Computing ',state,'(',irow,',',icol,')')
                # Getting the scatter data
                data = mxCOVID19.getFitData(state)
                # Getting the best fitted curve under the given hypothesis
                # Setting pars
                mxCOVID19.parVals = getPars(row)
                yValsModel = mxCOVID19.getModel(xtime)[:,1]
                # Analysis of the stat data
                mxCOVID19.getStatData(state)
                allPars = mxCOVID19.actStatData.copy()
                statData = allPars['eta'].values
                limitsETA = st.t.interval(conf, len(statData)-1, loc=np.mean(statData), scale=st.sem(statData))
                # Initializing the min and max y limits for the time range
                ymin = np.zeros(maxX)+1e8; ymax = np.zeros(maxX)
                minR0 = 1e5; maxR0 = 0; cnt = 0
                # Cycle over the statpars
                for i in range(len(allPars.index)):
                    statpars = allPars.iloc[i].values
                    if ((not useconf) or (statpars[2]>=limitsETA[0]) and (statpars[2]<=limitsETA[1])):
                        minR0 = minR0 if (minR0<statpars[3]) else statpars[3]
                        maxR0 = maxR0 if (maxR0>statpars[3]) else statpars[3]
                        mxCOVID19.parVals = getPars(statpars,True)
                        mxCOVID19.fE0 = statpars[9]
                        yVals = mxCOVID19.getModel(xtime)[:,1]
                        # Running over x to update max and min values
                        cnt+=1
                        for xindx in range(maxX):
                            ymin[xindx] = ymin[xindx] if (ymin[xindx]<yVals[xindx]) else yVals[xindx]
                            ymax[xindx] = ymax[xindx] if (ymax[xindx]>yVals[xindx]) else yVals[xindx]
                print(conf*100,'% on eta of '+state+':',limitsETA,'. Min R0=',minR0,', Max R0=',maxR0,' useconf=', useconf,' pars used: ',cnt)
                outres.append([state,np.max(ymin),xdate.iloc[np.argmax(ymin, axis=0)],minR0,maxR0])

        if (nstate<(len(source.index)+maxstates)):
            # The main axis
            # Plotting the best fit
            axs[irow,icol].plot(xdate,yValsModel,linewidth=0.8,color=colors[poscol])
            # Plotting the shadow probabilities
            axs[irow,icol].fill_between(xdate, ymin, ymax, alpha=0.3,color=colors[poscol])
            axs[irow,icol].tick_params(labelrotation=22)
            axs[irow,icol].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            upv = 0.05; starttext = 0.10; sz = 6; xpos = 0.67
            axs[irow,icol].text(xpos,starttext+upv*2, r'$\beta$='+'{:.2f}'.format(row[1]), transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(xpos,starttext+upv*1, r'$\eta$='+'{:.2f}'.format(row[2]), transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(xpos,starttext+upv*0, r'$R_0$='+'{:.1f}'.format(row[3]), transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(xpos,starttext+upv*4, r'$I_{m}$='+format(int(row[5]), ',d'), transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(xpos,starttext+upv*3, r'$D_{m}$='+row[6], transform=axs[irow,icol].transAxes,size = sz)
            axs[irow,icol].text(0.45,0.95, state, transform=axs[irow,icol].transAxes)#,color='blue', transform=ax.transAxes)

            # The inset axis
            axins.append(inset_axes(axs[irow,icol], width="30%", height="40%" ,borderpad=1.1))#, loc='lower left', bbox_to_anchor=(0.8, 0.8, 1, 1)))
            # Plotting the data
            axins[-1].scatter(mxCOVID19.realDays.values,data['Infected'],color='gray',s=2)
            # Plotting the shadow probabilities
            axins[-1].fill_between(xtime, ymin, ymax, alpha=0.3,color=colors[poscol])
            # Plotting the best fit
            axins[-1].plot(xtime,yValsModel,linewidth=0.8,color=colors[poscol])
            axins[-1].set_xlim(mxCOVID19.realDays.values[0],mxCOVID19.realDays.values[-1]) #self.fitData.iloc[0].values
            axins[-1].set_ylim(0,np.max(data['Infected'].values)) #self.fitData.iloc[0].values
            axins[-1].tick_params(direction='out', labelsize = 6)
            axins[-1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            print('Maximum for the most optimistic scenario: ', np.max(ymin), " achieved on ", xdate.iloc[np.argmax(ymin, axis=0)])

            # Moving to the next state
            nstate += 1

outres = pd.DataFrame(outres,columns = ['ID','MaxInf','Date','R0min','R0max'])
outres.to_csv (rawData[:8]+'_bestScenario.csv', index = False, header=True)
            
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.ylabel('Infected')
plt.tight_layout()
plt.savefig('FigureFITstates.png',dpi=600)
plt.show()

