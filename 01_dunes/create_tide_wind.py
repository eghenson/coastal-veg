import numpy as np

# this script creates two files:
# 
# a wind-file with three columns of values (wind.txt)
# the first column is a list of seconds, with a step of 3600 s (=1 h) for 31536000 s (=1 year)
# the second column is a list of wind speed, for now constant at 10 m/s
# the third column is a list of wind direction, for now constant at 0 degrees

# a tide-file with two columns of values (tide.txt)
# the first column is a list of seconds, with a step of 600 s (=10 min) for 31536000 s (=1 year)
# the second column is a list of tide height: a harmonic function for a semi-diurnal tide with a mean height of 1.5 m

# wind
t = np.arange(0, 31536000, 3600)
wind_speed = np.full_like(t, 10)
wind_dir = np.full_like(t, 0)
wind = np.column_stack((t, wind_speed, wind_dir))
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\wind.txt', wind, fmt='%i %f %f')

# tide
t = np.arange(0, 31536000, 600)
tide_height = 1.5 * np.sin(2 * np.pi * t / 44712)
tide = np.column_stack((t, tide_height))
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\tide.txt', tide, fmt='%i %f')

# Plotting
import matplotlib.pyplot as plt
plt.plot(t, tide_height)
plt.show()

