# This file only selects the middle three transects from the total 2D grid
# Input files:
# c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\x_old.grd
# c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\y_old.grd
# c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\z_old.grd
# c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\veg_old.grd
# c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\ne_old.grd

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Read in the x, y, z, veg, and ne grids
x_o = np.loadtxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\x_old.grd')
y_o = np.loadtxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\y_old.grd')
z_o = np.loadtxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\z_old.grd')
veg_o = np.loadtxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\veg_old.grd')
ne_o = np.loadtxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\ne_old.grd')

# Select the middle three transects
x = x_o[5, :]
y = y_o[5, :]
z = z_o[5, :]
veg = veg_o[5, :]
ne = z_o[5, :] - 0.011

# Write out the new grids
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\x.grd', x)
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\y.grd', y)
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\z.grd', z)
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\veg.grd', veg)
np.savetxt(r'c:\Users\weste_bt\Github\csdms-coastal-vegetation\01_building_dunes\setup\ne.grd', ne)


