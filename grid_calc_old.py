# Jarrad Pond
# SL9 Impact Modeling Team
# Unversity of Central Florida
# gridcalc.py script

#########
## This script, given minimal information, calculates the FORTRAN code
## Need to construct a computational grid following desired grid constraints
##
#########

# REAQUIRES (in same directory): read_atm.py, drake_jovian_atmo.txt

#####
## CODE START
#####

# import nessecary packages
import numpy as np

#####
## GIVE INPUTS
#####

# Diameter of impactor
D = 1000.e2 #cm

# Desired resolution across impactor diameter
R = 32 #elements per impactor diameter

# Calculate needed element/distance ratio needed for 
# uniform grid calculations
r = R/D # elements per cm of grid block

# Initial Absolute Height of impactor
#x10 = 2.315e07 #cm
z0 = 150.00 # km
z0 *= 1.e5 #cm

# Angle
angle = 45. #degrees from verical
angle_rad = angle * np.pi / 180.0 # change to radians
cosang = np.cos(angle_rad) # Take Cosine
sinang = np.sin(angle_rad) # Take Sine
print "Angle: %.6f" %angle_rad

#print angle_rad

# Calculate starting x1 height.
#z0 = x10 * cosang # cm
#print "z0 = %.3f km" %(z0/1e5)
x10 = z0 / cosang # cm
print "x10 = %.6e cm" %x10


# Read in atmosphere information
#from read_atm import *

# Find maximum absolute height of the atmosphere.
# z_atm_max = np.max(z_atm)
z_atm_max = 5.0e07 / 1e5 #km
print "Maximum atmosphere height: %.3f km \n" %(z_atm_max) 


# UNIFORM GRID CONTAINING PLUME INPUTS
# This grid can be constructed to be basically however you want it,
# and is mainly adjusted by giving the grid spacing and number of 
# elements in that space.

# N tile information
ntiles1 = 3.

ntiles2 = 2.

ntiles3 = 2.


# Grid spacing ratios.

# x1rat refers to the ratios for the x1 direction.
# there should be 4 here: 
# a ratio (non-1) for the buffer zone
# a ratio of 1.0 for the uniform impactor grid
# a ratio of 1.035 (or less than 4% change) for after the uniform impactor grid
# a raito of 1.0 for the uniform plume grid
x1rat = np.array([1.076, \
                  1.000, \
                  1.035])

# x2rat refers to the ratios for the x2 direction.
# there should be 3 here:
# a ratio (non-1) for the buffer zone on one side of the impactor
# a ratio of 1.0 for the uniform impactor grid
# a ratio (non-1) for the buffer zone on the other side of the impactor
x2rat = np.array([1.05, \
                  1.00, \
                  1.05])

# x3rat refers to the ratios for the x3 direction.
# there should be 3 here:
# a ratio (non-1) for the buffer zone on one side of the impactor
# a ratio of 1.0 for the uniform impactor grid
# a ratio (non-1) for the buffer zone on the other side of the impactor
# This is generally the same as x2rat
x3rat = np.array([1.05, \
                  1.00, \
                  1.05])


#####
## MAKE X1 GRID
#####

#####
## Calulate uniform impactor grid
#####

x12min = x10 - (1.5*D)
x12max = x10 + (2.5*D)

dx12 = x12max - x12min

nbl12 = np.int(np.round(r*dx12))

dx12min = dx12 / nbl12

#####
## Calulate non-uniform grid before impactor
#####

x11min = x10 - ((0.5)*D) - 9500.e2
x11max = x10 - (1.5*D)

dx11 = x11max - x11min

# We want dx11min to be within 3% of dx12min,
# so, we iterate ober nbl11 until we get a dx11min
# that is within 3% of dx12min

#ref = 1.0 - ((dx11*(1-x1rat[0]))/(dx12min*0.989))

x1rat0 =  x1rat[0]

x11mintest = x12min
for i in range(1,100):
	x11mintest = x11mintest - (dx12min*(x1rat0**(i-1)))
	if (x11mintest <= x11min):
		nbl11 = i
		x11min = x11mintest
		break
	
#####
## Calulate non-uniform grid after impactor
#####

x13min = x10 + (2.5*D)

# Now, we only want this grid to extend until dx13max is within
# 3% of dx14min. So, we will pick a dx13min that is within 3% of dx12min,
# increase x13 and dx13 until we get to a dx13 within 3% of dx14min.

# Initialize dx13min and x13
dx13min = dx12min*1.00
x13 = x13min

# Give total distance of grid
dx1_total = 3.6464e7 #cm

x13max = x11min + dx1_total
print(x13max)
dx13 = x13max - x13min

x1rat2 =  x1rat[2]

x13maxtest = x12max
for j in range(1,200):
	x13maxtest = x13maxtest + (dx13min*(x1rat2**(j-1)))
	if (x13maxtest >= x13max):
		nbl13 = j
		x13max = x13maxtest
		break


#####
## PRINT OUT FORTRAN CODE
#####
fmtstr = "%.6e"
fmtstr2 = "%.6e"
print ' $ggen1 nbl='+str(nbl11)+',x1min='+str(fmtstr %x11min)+',x1max='+str(fmtstr %x11max)+',igrid=-1,x1rat='+str("%.4f" %x1rat[0])+',lgrid=.false. /'
print ' $ggen1 nbl='+str(nbl12)+',x1min='+str(fmtstr %x12min)+',x1max='+str(fmtstr %x12max)+',igrid=1,x1rat='+str("%.4f" %x1rat[1])+',lgrid=.false. /'
print ' $ggen1 nbl='+str(nbl13)+',x1min='+str(fmtstr %x13min)+',x1max='+str(fmtstr %x13max)+',igrid=1,x1rat='+str("%.4f" %x1rat[2])+',lgrid=.true. /'


#####
## X1 GRID MADE
#####


#####
## MAKE X2 AND X3 GRIDS
#####

# In most cases, these grid will be identical.

#####
## Calulate uniform impactor grid
#####

x22min = 0.0 - (1.0*D)
x22max = 0.0 + (1.0*D)

dx22 = x22max - x22min

nbl22 = np.int(np.round(r*dx22))

dx22min = dx22 / nbl22

#####
## Calulate non-uniform grid around impactor, in the negative direction.
#####

# x21min = 0.0 - ((1.0 + 101)*D)
# x21max = 0.0 - (1.0*D)

x21min = 0.0 - 100.e5 
x21max = 0.0 - (1.0*D)


dx21 = x21max - x21min

# We want dx21min to be within 3% of dx22min,
# so, we iterate ober nbl21 until we get a dx21min
# that is within 3% of dx22min

#ref = 1.0 - ((dx11*(1-x1rat[0]))/(dx12min*0.989))

x2rat0 = x2rat[0]

x21mintest = x22min
for k in range(1,120):
	x21mintest = x21mintest - (dx22min*(x2rat0**(k-1)))
	if (x21mintest <= x21min):
		nbl21 = k
		x21min = x21mintest
		break

#####
## Calulate non-uniform grid around impactor, in the positive direction.
#####
x23min = -1.0*x21max
x23max = -1.0*x21min
nbl23 = nbl21

#####
## PRINT OUT FORTRAN CODE
#####
# x2
print ' $ggen2 nbl='+str(nbl21)+',x2min='+str(fmtstr %x21min)+',x2max='+str(fmtstr %x21max)+',igrid=-1,x2rat='+str("%.4f" %x2rat[0])+',lgrid=.false. /'
print ' $ggen2 nbl='+str(nbl22)+',x2min='+str(fmtstr %x22min)+',x2max='+str(fmtstr %x22max)+',igrid=1,x2rat='+str("%.4f" %x2rat[1])+',lgrid=.false. /'
print ' $ggen2 nbl='+str(nbl23)+',x2min='+str(fmtstr %x23min)+',x2max='+str(fmtstr %x23max)+',igrid=1,x2rat='+str("%.4f" %x2rat[2])+',lgrid=.true. /'
# x3 is the same
print ' $ggen3 nbl='+str(nbl21)+',x3min='+str(fmtstr %x21min)+',x3max='+str(fmtstr %x21max)+',igrid=-1,x3rat='+str("%.4f" %x2rat[0])+',lgrid=.false. /'
print ' $ggen3 nbl='+str(nbl22)+',x3min='+str(fmtstr %x22min)+',x3max='+str(fmtstr %x22max)+',igrid=1,x3rat='+str("%.4f" %x2rat[1])+',lgrid=.false. /'
print ' $ggen3 nbl='+str(nbl23)+',x3min='+str(fmtstr %x23min)+',x3max='+str(fmtstr %x23max)+',igrid=1,x3rat='+str("%.4f" %x2rat[2])+',lgrid=.true. /'


#######
## GRID BACK CHECK
## Check if the back of the grid, in the initial configuration, 
## sticks out of the top of the simulated atmosphere.
#####

# Caculate the absolute height of the top corner of the grid.
z_top_corner = (x13max * cosang) + (x23max * sinang) #cm
z_top_corner /= 1e5 #km
print "\nTop corner location = %.3f km" %z_top_corner

if (z_top_corner > z_atm_max):
    print "Grid Bad. Lower starting height"
    print "Difference = ", (z_top_corner - z_atm_max)

if (z_top_corner < z_atm_max):
    print "Grid Fine!"
    print "Difference = ", (z_top_corner - z_atm_max)

#####
## DX1SAFE CALCULATION
#####

# Calculate the appropriate dx1safe (the distance between
# the front of the impactor and the bottom of the grid)
# so that impactor doesn't fall out of the grid.

dx1safe = x10 - D/2.0 - x11min
print 'Use dx1safe = %.3e' %dx1safe
print 'Grid Length: dx1 = %.5e' %(x13max-x11min) + ' dx2 =  %.5e' %(x23max-x21min) + ' dx3 =  %.5e' %(x23max-x21min) 

print 'r = %.3e cm' %(D/2.0)

#####
## NODE CALCULATIONS
#####

print "\n NODE CALCULATIONS:"

# Print # or nodes used
print "Total Number of nodes used: " + str(ntiles1) + " x " + str(ntiles2) + " x " +str(ntiles3) + " = " + str((ntiles1 * ntiles2 * ntiles3))

# Print the total number of cells in each direction
print "Number of cells in x1: ", (nbl11 + nbl12 + nbl13)
print "IZONES = ", ((nbl11 + nbl12 + nbl13)/ntiles1)
print "Number of cells in x2: ", (nbl21 + nbl22 + nbl23)
print "JZONES = ", ((nbl21 + nbl22 + nbl23)/ntiles2)
print "Number of cells in x3: ", (nbl21 + nbl22 + nbl23)
print "KZONES = ", ((nbl21 + nbl22 + nbl23)/ntiles3)
