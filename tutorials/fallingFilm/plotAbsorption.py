### IMPORTS

import sys
import re
import matplotlib.pyplot as plt
import math
import numpy as np

### MAIN

# Set basic plot parameters
fig, (ax1) = plt.subplots(1, sharex=False)

# Font size
fs = 16

ax1.set_ylabel("CO2 Absorption (kg/s)", fontsize=fs)
ax1.set_xlabel("Time (s)", fontsize=fs)

for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(14)

# Post-process each log file given in command line argument
for q in range(len(sys.argv)-1) :

    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()
    
    # Initialize empty data arrays
    mDotBulkIn = []
    mDotBulkOut = []
    mDotFilmIn = []
    mDotFilmOut = []
    times = []
    
    # Search log file line-by-line for desired data
    for line in lines :
        	
        if "Bulk mDotCO2in" in line :
            mDotBulkIn.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[1]))
            
        if "Bulk mDotCO2out" in line :
            mDotBulkOut.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[1]))
            
        if "Film mDotCO2in" in line :
            mDotFilmIn.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[1]))
            
        if "Film mDotCO2out" in line :
            mDotFilmOut.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[1]))
            
        if "Time" in line and not "Execution" in line :
            times.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[0]))
            
    times.pop(0)
    mDotBulkIn.pop(0)
    mDotBulkOut.pop(0)
    mDotFilmIn.pop(0)
    mDotFilmOut.pop(0)
    
    # Get net change
    mDotBulkIn = np.array(mDotBulkIn)
    mDotBulkOut = np.array(mDotBulkOut)
    mDotBulk = -np.add(mDotBulkIn, mDotBulkOut)
    mDotFilmIn = np.array(mDotFilmIn)
    mDotFilmOut = np.array(mDotFilmOut)
    mDotFilm = np.add(mDotFilmIn, mDotFilmOut)
    
    if (len(times) > len(mDotFilm)) :
        times.pop()

    # Plot changes
    ax1.plot(times, mDotBulk, lw=2.5, label="Bulk")
    ax1.plot(times, mDotFilm, lw=2.5, label="Film")
    
# Plot target bulk velocity
#ax1.axhline(y = 3.26, color = 'black', linestyle = '-.', lw=2.5, label="Target Bulk Velocity")

# Set miscellaneous plot parameters
ax1.legend(prop={'size': 14})
ax1.grid(True)
#ax1.set_xlim(xmin=0)
#ax1.set_ylim(ymin=0)

plt.subplots_adjust(left=0.15,
                    bottom=0.12,
                    right=0.95,
                    top=0.95,
                    wspace=0.2,
                    hspace=0.4)

# Show plot window
plt.show()
