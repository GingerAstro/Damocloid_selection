# importing apps
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
pd.set_option('display.max_columns', 500)

print('Reading in data')
# read in MPCORB.DAT
col = [(0,7),(8,13),(14,19), (20,25), (26,35), (37,46), (48,57), (59,68), (70,79), (80,91),(92,103), (105,106),(107,116), (117,122), (123,126), (127,131), (132,136), (137,141), (142,145), (146,149), (150,160), (161,165), (166,194), (194,202)]
names = ['Number','h','g', 'epoch', 'anomaly', 'arg peri', 'node', 'i', 'e', 'daily motion', 'a', 'u', 'ref', 'num obs', 'num oppos', 'year first', 'year last', 'rms', 'coarse', 'precise', 'computer', 'hex flag', 'read desig', 'last obs']
asteroid = pd.read_fwf ('MPCORB.DAT', colspecs=col, names=names, header=None, skiprows=43)

# read in CometEls.txt
col = [(0,4),(4,5),(5,12), (14,19), (20,22), (22,30), (31,40), (41,50), (51,60), (61,70),(71,79), (81,85), (85,87), (87,89), (91,95), (96,100), (102,158), (159,168)]
names = ['Number','orbit','designation', 'year', 'month', 'day', 'peri', 'e', 'arg peri', 'node', 'i', 'year ep', 'month ep', 'day ep', 'h', 'slope', 'name', 'ref']
comet = pd.read_fwf ('CometEls.txt', colspecs=col, names=names, header=None, skiprows=0)

print('Data read finished')
print('Calculating paramaters')
# filter out incomplete data
comet = comet[np.isfinite(comet['peri'])]
comet = comet[np.isfinite(comet['e'])]
comet = comet[comet.e<1] # removing hyperbolic comets

asteroid = asteroid[np.isfinite(asteroid['a'])]
asteroid = asteroid[np.isfinite(asteroid['e'])]

# Period calculation
comet['a'] = comet['peri']/(1-comet['e'])
comet['period'] = np.sqrt(comet['a']**3)
comet['aph'] = comet['a']*(1+comet['e'])

asteroid['period'] = np.sqrt(asteroid['a']**3)
asteroid['aph'] = asteroid['a']*(1+asteroid['e'])
asteroid['peri'] = asteroid['a']*(1-asteroid['e'])

# Tisserand parameter calculations
aj= 5.204

comet['tiss'] = (aj/comet['a'])+(2*(np.sqrt((comet['a']/aj)*(1.0-comet['e']**2)))*np.cos(np.radians(comet['i'])))

asteroid['tiss'] = (aj/asteroid['a'])+(2*(np.sqrt((asteroid['a']/aj)*(1.0-asteroid['e']**2)))*np.cos(np.radians(asteroid['i'])))

print('Calculations compleated')

# seperation of long and short period comets
lpcomet = comet[comet.period>200]
spcomet = comet[comet.period<200]

# calculate damocloids in different methods
jewitt = asteroid[asteroid.tiss<2]
gibson = asteroid[(asteroid.aph>65) & (asteroid.e>0.9)]

#saving damocloid targets
datestring = datetime.strftime(datetime.now(), ' %Y_%m_%d')
gibson.to_csv('gibson_damocloids' + datestring + '.csv')

print('Number of Gibson damocloids detected = ', len(gibson))
print('Number of Jewitt damocloids detected = ', len(jewitt))
print('Total Comets = ', len(comet))
print('Total lpComets = ', len(lpcomet))
print('Total spComets = ', len(spcomet))