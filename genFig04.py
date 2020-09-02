import geopandas as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Reading mexico map shape
states = gp.geopandas.read_file("mex/mexican-states.shp")
res20200525 = pd.read_csv("20200525_ResXStates.csv")

# Dictionary to convert three leter codes for states names (some of the codes included in the dataframe are different)
statesDict = {'Aguascalientes':'AGU','Colima':'COL','Tlaxcala':'TLA','Ciudad de México':'CMX','Morelos':'MOR','México':'MEX','Hidalgo':'HID','Puebla':'PUE','Nuevo León':'NLE','Coahuila de Zaragoza':'COA','Chihuahua':'CHH','Sonora':'SON','Michoacán de Ocampo':'MIC','Querétaro':'QUE','Guanajuato':'GUA','Jalisco':'JAL','Zacatecas':'ZAC','Durango':'DUR','Tamaulipas':'TAM','Veracruz de Ignacio de la Llave':'VER','Guerrero':'GRO','Sinaloa':'SIN','Oaxaca':'OAX','Nayarit':'NAY','Chiapas':'CHP','Tabasco':'TAB','Campeche':'CAM','Baja California':'BCN','Baja California Sur':'BCS','San Luis Potosí':'SLP','Yucatán':'YUC','Quintana Roo':'ROO'}

# Getting R0 values
states['R0'] = [ res20200525.loc[res20200525['ID'] == statesDict[item]].values[0][3] for item in states.name.values]
# Setting plotting style
plt.style.use('seaborn')

################################################################################################
# R0 map plot
################################################################################################

vmin = states['R0'].min(); vmax = states['R0'].max();
states.plot(column='R0', cmap='Blues', norm=plt.Normalize(vmin=vmin,vmax=vmax), linewidth=0.3, edgecolor='.3')

sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
plt.colorbar(sm)

plt.title(r'$R_0$ values in Mexico', fontdict={'fontsize': '18', 'fontweight' : '3'})
plt.annotate('Reproductive number obtained from model fit up to 25/06/2020',xy=(0.2, .15),  xycoords='figure fraction', horizontalalignment='left', verticalalignment='top', fontsize=10, color='#555555')
plt.axis('off')
plt.savefig('Figure04.png',dpi=600)
plt.show()

