# Noah Nodolski and Ramanakumar Sankar 5/2/2018
import numpy as np
import sys




############## x1 gridgen ################




#### PARAMETERS #####################################
#       ncpus - number of cpus                                  #
#       s - size of impactor (cm)                               #
#       sres - desired resolution at the impactor               #
#       lbef - desired length of grid before the                #
#                       impactor with same resolution           #
#                       as impactor                             #
#       laft - desired length of grid after the                 #
#                       impactor with same resolution           #
#       l1orig - desired length of grid before the              #
#                       impactor with increasing grid spacing   #
#       l2orig - desired length of grid after the               #
#                       impactor with increasing grid spacing   #
#       maxr - maximum resolution of grid                       #
#       start - grid starting point                             #
#       All lengths in centimeters                              #
#                                                               #
#####################################################

if(sys.version_info < (3,0)):
        print("Correcting for Python 2.7")
        def input(str=""):
                return raw_input(str)

def get_inp(default):
        inp = input()
        if(inp == ""):
                inp = default
        return float(inp)

def get_inp_tuple(default1, default2):
        tup = input()
        if(tup == ""):
            tup = (default1, default2)
        else:
            tup = tuple(map(float,tup.split(' ')))
        return (tup)

def get_nbl(xmin, xmax, dx, r):
    dlength = xmax - xmin
    N = np.log((dlength*(r-1.) + dx)/dx)/\
        np.log(r)
    return int(np.ceil(N))


print("\n\n\n")
print('\n=============================================')
print(' Initial conditions and x1 grid generation   ')
print('=============================================\n')

good = False

while(good == False):
        ## Define maximum threshold with nzones and ncpus (only for com parison)
        # nzones = 62 ## UNNECESSARY
        print("Number of CPUs in this direction [Default 6]: ")
        ncpu = get_inp(6)
        # ntot = nzones*ncpu ## UNNECESSARY

        ## Define maximum resolution
        print("Maximum resolution (m per block) [4000]: ")
        maxr = get_inp(4000.)*1.e2

        ## Impactor size
        print("Size of impactor (m) [1000]: ")
        s = get_inp(1000.)*1.e2

        ## Impactor Velocity
        print("Velocity of Impactor (km/s) [61.4]")
        v = get_inp(61.4)*1.e5

        ## Impactor Viscosity
        print("Viscosity of Impactor [0.5]")
        qcon = get_inp(0.5)

        ## Resolution at the impactor
        sres = (s/32)

        ## Number of blocks target takes
        ntarget = np.floor(s/sres)

        ## Calculate the refined resolution
        sres = s/ntarget

        ## Define the length before and after target with same res as target
        lbef = 1.*s
        laft = 2.*s

        ## Starting point
        print("Starting height(km) [150]: ")
        h = get_inp(150.)*1.e3*1.e2

        print("Impact angle(deg) [45]: ")
        ang = get_inp(45.)*np.pi/180.

        ## Add the estimated grid size before and after l1 - before and l2 - after
        print("Estimated grid size before and after high res zone (m) [9500 and 364600]")
        lorig = get_inp_tuple(9500, 364600)

        l1orig = lorig[0]*1.e2
        l2orig = lorig[1]*1.e2


        start = h/np.cos(ang)

        ## Calculate number of blocks and length of grid before target
        nbef = np.floor(lbef*(ntarget/s))
        lbef = nbef*s/ntarget

        ## Do the same for after the target
        naft = np.floor(laft*(ntarget/s))
        laft = naft*s/ntarget

        ## Number of blocks for collider
        ncollider = ntarget + nbef + naft

        print("Collider zone: %d blocks, from %.4f m before to %.4f m after"%(ncollider,lbef/100., laft/100.))
        #print("Total block remaining: %d"%(ntot - ncollider))

        print("\n")

        ## Add to the blocks array
        blocks = [[nbef + ntarget + naft, 1., 1., lbef + s + laft]]

        blocks2 = []
        ## Read in ratios
        print("Ratio before, after [Default 1.05 and 1.03]: ")
        r = get_inp_tuple(1.05, 1.03)

        print("\t #Blks\t Ratio \t \t Length (m) \t Given length (m)")
        print("---------------------------------------------")

        ## get the number of blocks based on the initial high res
        ## grid size
        n1 = get_nbl(-l1orig, 0., sres, r[0])
        l1 = sres*(1.-(r[0]**(n1)))/(1.-r[0])
        l1best = l1
        scheck1 = l1*(r[0]-1.)/(r[0]**n1 - 1.)

        if(np.abs(sres-scheck1)/sres > 0.03):
                print("resolution mismatch before impactor ---- check numbers")

        n2 = get_nbl(0., l2orig, sres, r[1])
        l2 = sres*(1.-(r[1]**(n2)))/(1.-r[1])

        ## Calculate the maximum grid spacing at the endpoints
        last1 = sres*r[0]**(n1)
        last2 = sres*r[1]**(n2)

        total = 0.
        ''' FIGURE OUT LOW RES n BEFORE '''
        if(last1 > maxr):       ## If grid before has lower resolution than max
                ## Calculate the number of blocks to reach max resolution
                nmax1 = np.floor(np.log(maxr/sres)/np.log(r[0]))

                ## Calculate grid size upto block with max resolution
                lmax1 = sres*(r[0]**(nmax1-1) - 1.)/(r[0] - 1.)

                ## Calculate the max resolution
                rmax1 = sres*r[0]**(nmax1-1)

                blocks2.append([nmax1, -1., rmax1, lmax1])

                ## Calculate the number of flat (1.0 ratio) blocks to
                ## fill the remaining grid spacing
                nflat1 = np.floor((l1orig - lmax1)/rmax1)

                total = total + nflat1 + nmax1

                ## Calculate grid size of flat section
                lflat1 = rmax1*nflat1
                blocks2.append([nflat1,1., 1., lflat1])
                print("Before:\t %d \t %.4f \t %.4f"%(nmax1, r[0], lmax1/100.))
                print("\t\t %d \t %.4f \t %.4f"%(nflat1, 1., lflat1/100.))
                print("------------------------------------------")
                print("\t\t %d \t \t %.4f \t %.4f"%(nflat1 + nmax1, (lflat1 + lmax1)/100., l1orig/100.))
        else:   ## If max resolution is not exceeded
                blocks2.append([n1,-1., r[0], l1])
                total = total + n1
                print("Before:\t %d \t %.4f \t %.4f \t %.4f"%(n1, r[0],l1/100., l1orig/100.))

        ## Add the impactor grid
        blocks2.append(blocks[0])

        ## Rinse and repeat for the grid behind impactor
        print("\n")
        if(last2 > maxr):
                ## search for max r
                ltry = 0.
                for i in range(200):
                        ## there is a factor of 1.03 here... not sure why
                        rtry = 1.03*sres*(r[1]**float(i))
                        ltry += rtry
                        percdiff = (rtry - maxr)/maxr
                        if((percdiff<=0.03)&(percdiff>=0.)):
                                ## if the %diff is good then ouput the previous value
                                ## -1 because of this overestimates it
                                nmax2 = i - 1
                                lmax2 = ltry - rtry
                                rmax2 = maxr

                print("Cutting zone after at %d - max resolution %.3e"%(nmax2, rmax2))

                nmax2 = get_nbl(0., lmax2, sres, r[1])
                lmax2 = sres*(1.-(r[1]**(nmax2)))/(1.-r[1])
                rmax2 = sres*(r[1])**(nmax2-1)

                blocks2.append([nmax2, 1., r[1], lmax2])

                ## find the number of flat zones
                nflat2 = np.floor((l2orig - lmax2)/rmax2)
                nzones = (total + nmax2 + nflat2 + ncollider)/ncpu

                while(nzones%1. != 0.):
                        nflat2 = nflat2 + 1
                        nzones = (total + nmax2 + nflat2 + ncollider)/ncpu

                ## calculate the flat area
                lflat2 = rmax2*(nflat2)
                blocks2.append([nflat2, 1., 1., lflat2])
                total = total + nflat2 + nmax2
                print("After:\t %d \t %.4f \t %.4f"%(nmax2, r[1], lmax2/100.))
                print("\t %d \t %.4f \t %.4f"%(nflat2, 1., lflat2/100.))
                print("------------------------------------------")
                print("\t %d \t \t \t %.4f \t %.4f"%(nflat2 + nmax2, (lflat2 + lmax2)/100., l2orig/100.))
                l2best = lmax2 + lflat2
        else:
                nzones = (total + n2 + ncollider)/ncpu

                while(nzones%1. != 0.):
                        n2 = n2 + 1
                        nzones = (total + n2 + ncollider)/ncpu

                while(nzones%1. != 0.):
                        n2 = n2 + 1
                        nzones = (total + n2 + ncollider)/ncpu

                l2 = sres*(1.-(r[1]**(n2)))/(1.-r[1])
                l2best = l2

                blocks2.append([n2, 1., r[1], l2])
                total = total + n2
                print("After:\t %d \t %.4f \t %.4f \t %.4f"%(n2, r[1],l2/100., l2orig/100.))

        ## Print out total number of blocks and number of zones (from ncpus)
        print("Total:\t %d "%(total))
        print("Total blocks: %.2f in %.2f zones"%(total+ncollider,(total+ncollider)/ncpu))

        ## calculate the block of the front and back of impactor
        ## to check where zone division takes place
        impfronti = blocks2[0][0] + 32
        impbacki  = blocks2[0][0] + 64

        izonetest = (total+ncollider)/ncpu

        zonebef = 0
        while(zonebef < impfronti):
                zonebef += izonetest
        zonebef = zonebef - izonetest

        print("Impactor front, Zone division: %d %d"%(impfronti, zonebef))

        zoneaft = zonebef
        while(zoneaft <= impbacki):
                zoneaft += izonetest

        print("Impactor back, Zone divison: %d %d"%(impbacki, zoneaft))

        ans = input("Are these values correct? (y/n): ")

        if(ans in ['y','Y','yea','yes','Yes','Yeah','yeah']):
                good = True
                blocks = blocks2
        else:
                good = False

izones = (total+ncollider)/ncpu




############## x2 and x3 gridgen ################



print('\n===============================================')
print('            x2 and 3 grid generation           ')
print('===============================================\n')

# set loop condition

bad = True
jzones = 3.5
kzones = 3.5

# Calculate Low-res distances
def zones(d,r,b):
        length = (d/32)*((1.-r**b)/(1.-r)) + d
        return length

while(bad == True):

        jzones = 3.5
        kzones = 3.5

        while((jzones %1 != 0) or (kzones %1 != 0)):
                # Define number of cpus for zone calculations
                print("Enter number of cpus for x2 and x3 direction [Default 2 and 2]: ")
                ncpu2 = get_inp_tuple(2, 2)
                
                print("Enter the grid length before and after impactor for x2 dimension in km [Default %.4f km and %.4f km]: " %((l2/100000)*0.4, (l2/100000)*0.4))
                x2_len = get_inp_tuple(l2*0.4, l2*0.4)
                if ((x2_len[0] != l2*0.4) & (x2_len[1] != l2*0.4)):
                    x2_len = (x2_len[0]*1.e5, x2_len[1]*1.e5)

                # Length Before and after Impactor
                print("Enter the grid length before and after impactor for x3 dimension in m [Default %.4f km and %.4f km]: " %((l2/100000)*0.4, (l2/100000)*0.4))
                x3_len = get_inp_tuple(l2*0.4, l2*0.4)
                if ((x3_len[0] != l2*0.4) & (x3_len[1] != l2*0.4)):
                    x3_len = (x3_len[0]*1.e5, x3_len[1]*1.e5)

                # Ratios
                print("Ratio before, after for x2 [Default 1.05 and 1.05]: ")
                x2_r = get_inp_tuple(1.05, 1.05)

                print("Ratio before, after for x3 [Default 1.05 and 1.05]: ")
                x3_r = get_inp_tuple(1.05, 1.05)

                x2_blk = np.zeros(2)
                x3_blk = np.zeros(2)
                x2_blk[0] = get_nbl(-x2_len[0], -s, sres, x2_r[0])
                x2_blk[1] = get_nbl(s, x2_len[1], sres, x2_r[1])

                x3_blk[0] = get_nbl(-x3_len[0], -s, sres, x3_r[0])
                x3_blk[1] = get_nbl(s, x3_len[1], sres, x3_r[1])

                # Calculate Low-res zones before and after impactor in the x2 direction
                x2_min1 = zones(s, x2_r[0], x2_blk[0]) * -1
                x2_max3 = zones(s, x2_r[1], x2_blk[1])

                # Calculate Low-res zones before and after impactor in the x3 direction
                x3_min1 = zones(s, x3_r[0], x3_blk[0]) * -1
                x3_max3 = zones(s, x3_r[1], x3_blk[1])

                # Calculate total number of blocks and zones of each dimension
                jblocks = x2_blk[1] + 64 + x2_blk[0]
                kblocks = x3_blk[1] + 64 + x3_blk[0]
                jzones = (jblocks)/ncpu2[0] 
                kzones = (kblocks)/ncpu2[1]
                
                if((jzones %1 != 0) or (kzones %1 != 0)):
                        print('\n\n********************************************************')
                        print('Non integer number of zones. Please Re-enter Parameters.') 
                        print('********************************************************\n\n')

        # Output Calculations for looking over
        print("\n")
        print(("First x2 low-res zone: %.4f km\n") %(x2_min1*1.e-5*-1))
        print(("Second x2 low-res zone: %.4f km\n") %( x2_max3*1.e-5))
        print(("Total x2 blocks %.2f in %.3f zones\n\n") %(jblocks, jzones))
        print(("First x3 low-res zone: %.4f km\n") %(x3_min1*1.e-5*-1))
        print(("Second x3 low-res zone: %.4f km\n") %(x3_max3*1.e-5))
        print(("Total x3 blocks %.2f in %.3f zones\n\n") %(kblocks, kzones))

        ans = input("Are these values correct? (y/n):")
        if(ans in ['y','Y','yea','yes','Yes','Yeah','yeah']):
                bad = False
        else:
                bad = True

# define outputs frequency
print("\nEnter output frequency [Default 0.25]:")
dthdf = get_inp(0.25)

#define timeslice frequency
print("Enter timeslice frequency [Default 0.25]:")
dttsl = get_inp(0.25)



############ Write sim_log file ################ 

sim_name = input("\nEnter Simulation Name: ")

g = open('sim_log', 'w')
g.write("*********************************************************************\n")
g.write("\t\t\tImpact simulation %s\n"%(sim_name))
g.write("*********************************************************************\n\n")
g.write('Impactor Diameter = %.1f m\n\n'%(s/100))
g.write('Starting Velocity  = %.1f km/s\n\n'%(v/(1.e5)))
g.write('Viscosity = %.1f\n\n'%(qcon))
g.write('Starting Height = %.1f km\n\n'%(h/(1.e5)))
g.write('Entry Angle = %.1f degrees\n\n'%(ang*180/np.pi))
g.write('Length before impactor in x1 direction = %.1f km\n\n'%(l1best/(1.e5)))
g.write('Length after impactor in x1 direction = %.1f km\n\n'%(l2best/(1.e5)))
g.write('Length below impactor in x2 direction = %.1f km\n\n'%(x2_min1*1.e-5*-1))
g.write('Length above impactor in x2 direction = %.1f km\n\n'%(x2_max3*1.e-5))
g.write('Length below impactor in x3 direction = %.1f km\n\n'%(x3_min1*1.e-5*-1))
g.write('Length above impactor in x3 direction = %.1f km\n\n'%(x3_max3*1.e-5))
g.write('Processors configuration [(x1) x (x2) x (x3)] = %d x %d x %d\n\n'%(ncpu, ncpu2[0], ncpu2[1]))
g.write('Output frequency = %.2f\n\n'%(dthdf))
g.write('Timeslice frequency = %.2f\n\n\n'%(dttsl))
g.write('Comments:\n')

############ Write zmp_inp file ################ 




f = open('zmp_inp', 'w')

f.write(' &GEOMCONF  LGEOM    = 1,\n')
f.write('            LDIMEN   = 3 /\n')
f.write(' &PHYSCONF  LRAD     = 0,\n')
f.write('            LEOS     = 1,\n')
f.write('            NSPEC    = 3,\n')
f.write('            XHYDRO   = .TRUE.,\n')
f.write('            XFORCE   = .TRUE.,\n')
f.write('            XMHD     = .FALSE.,\n')
f.write('            XTOTNRG  = .false.,\n')
f.write('            XGRAV    = .false.,\n')
f.write('            XPTMASS  = .false.,\n')
f.write('            XISO     = .false.,\n')
f.write('            XSUBAV   = .false.,\n')
f.write('            XVGRID   = .true.,\n')
f.write('            XZGRAV   = .true. /\n')
f.write(' &IOCONF    XASCII   = .false.,\n')
f.write('            XHDF     = .true.,\n')
f.write('            XTSL     = .true.,\n')
f.write('            XFLX     = .true.,\n')
f.write('            XRESTART = .true.,\n')
f.write('            XPCLE    = .false. /\n')
f.write(' &PRECONF   SMALL_NO = 1.0D-99,\n')
f.write('            LARGE_NO = 1.0D+99 /\n')
f.write(' &ARRAYCONF IZONES   = %d ,\n'%(izones))
f.write('            JZONES   = %d,\n'%(jzones))
f.write('            KZONES   = %d,\n'%(kzones))
f.write('            MAXIJK   = %d /\n'%(np.max([izones,jzones,kzones])))
f.write(' &ATMCONF   XGND     = .false.,\n')
f.write('            XJUP     = .false.,\n')
f.write('            XJUP_DRAKE = .true.,\n')
f.write('            XTITAN   = .false.,\n')
f.write('            XVENUS   = .false. /\n')
f.write(' &mpitop ntiles(1)=%d,ntiles(2)=%d,ntiles(3)=%d,periodic=3*.false. /\n'%(ncpu, ncpu2[0], ncpu2[1]))
f.write(" &rescon tdump=0.0, dtdump=-0.5,id='ae',irestart=0,resfile='resae000000.000'/\n")
f.write(' &pcon nlim=9999999, tlim=20.0, cpulim= 9999999.0, mbatch=1 /\n')
f.write(' &hycon qcon=%.1f,courno=0.65,dtrat=1.0e-03 /\n'%(qcon))
f.write(' &iib niis(1)=6 /\n')
f.write(' &oib nois(1)=2 /\n')
f.write(' &ijb nijs(1)=2 /\n')
f.write(' &ojb nojs(1)=2 /\n')
f.write(' &ikb niks(1)=2 /\n')
f.write(' &okb noks(1)=2 /\n')



## Define grid starting point
offset = h/np.cos(ang) - blocks[0][3] - lbef - 0.5*s
# offset = h/np.cos(ang) - 1.5*s - lbef

total = 0.
block = []
done = 'false'

## Create the blocks text to output
for i in range(len(blocks)):
        n, dir, r, l = blocks[i]
        total = total + n
        if(i == (len(blocks) - 1)):
                done = 'true'
        x1 = offset
        x2 = offset + l

        if(i==0):
                startingpoint = x1
        elif(i==1):
                starthighres = x1

        offset = x2
        f.write("$ggen1 nbl=%d,x1min=%.4e,x1max=%.4e,igrid=%d,x1rat=%.3f,lgrid=.%s. /\n"%(n, x1, x2, dir, r, done))

dx1 = (starthighres + s) - startingpoint
x10 = starthighres + (3./2.)*s

f.write("$ggen2 nbl=%d,x2min=%.4e,x2max=%.4e,igrid=%d,x2rat=%.3f,lgrid=.false. /\n" %(x2_blk[0], x2_min1, s * -1, -1, x2_r[0]))
f.write("$ggen2 nbl=%d,x2min=%.4e,x2max=%.4e,igrid=%d,x2rat=%.3f,lgrid=.false. /\n" %(64, s * -1, s, 1, 1))
f.write("$ggen2 nbl=%d,x2min=%.4e,x2max=%.4e,igrid=%d,x2rat=%.3f,lgrid=.true. /\n" %(x2_blk[1], s, x2_max3, 1, x2_r[1]))
f.write("$ggen3 nbl=%d,x3min=%.4e,x3max=%.4e,igrid=%d,x3rat=%.3f,lgrid=.false. /\n" %(x3_blk[0], x3_min1, s * -1, -1, x3_r[0]))
f.write("$ggen3 nbl=%d,x3min=%.4e,x3max=%.4e,igrid=%d,x3rat=%.3f,lgrid=.false. /\n" %(64, s * -1, s, 1, 1))
f.write("$ggen3 nbl=%d,x3min=%.4e,x3max=%.4e,igrid=%d,x3rat=%.3f,lgrid=.true. /\n" %(x3_blk[1], s, x3_max3, 1, x3_r[1]))
f.write(' &grvcon /\n')
f.write(' &radcon /\n')
f.write(' &eqos gamma=1.4,mmw=1.0D0 /\n')
f.write(' &pgen r=%.4e, x10=%.4e , d0=0.917, e0=0.6, v11=%.2e, e1=0.0, \n'%(s/2., x10, v * -1))
f.write('     dx1safe=%.4e ,zpost0=-2.50e07,dzpost=4.0e07,theta0=%.4f,\n'%(dx1, ang))
f.write('     lat0=44.02 /\n')
f.write(' &gcon igcon=5,x1fac=1.0/\n')
f.write(' &iocon thdf=0.0, dthdf=%.2f, ttsl=0.0, dttsl=%.2f, tpcle=0.0, dtpcle=0.5 /\n'%(dthdf, dttsl))


f.close()


print ('\n\n***DONE***\n\n')


