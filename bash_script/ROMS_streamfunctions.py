import numpy
from netCDF4 import Dataset
import shutil
import argparse

from computeOSF import computeOSF

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--plot", dest="plot", action='store_true',
                    help="Perform plotting at each time step.")
parser.add_argument("-r", "--roms_file", dest="roms_file",
                    help="The ROMS input file", required=True)
parser.add_argument("-i", "--input_file", dest="input_file",
                    help="The ISOMIP+ input file", required=True)
parser.add_argument("-o", "--output_file", dest="output_file",
                    help="The ISOMIP+ output file", required=True)
args = parser.parse_args()

# first, make a copy of the input file in which we'll just overwrite the OSF
shutil.copy(args.input_file, args.output_file)

x0 = 320e3
y0 = 0.

dx = 2e3
dy = 2e3

x0ROMS = x0 + dx  # because ROMS has a 1 grid cell buffer around it
y0ROMS = y0 + dy  # because ROMS has a 1 grid cell buffer around it

nxOut = 240
nyOut = 40
nzOut = 144
dzOut = 720./nzOut

nzExtra = 1

# the z location of grid-cell corners (or top-bottom inderfaces) on the output
# grid, with an extra point at the top to catch anything above 0
zInterfaceOut = -dzOut*numpy.arange(-nzExtra, nzOut+1)
# the z location of grid-cell centers on the output grid
zOut = 0.5*(zInterfaceOut[0:-1] + zInterfaceOut[1:])

ncROMS = Dataset(args.roms_file, 'r')

variables = ncROMS.variables

nTime = len(ncROMS.dimensions['ocean_time'])

hc = variables['hc'][:]

s_rho = variables['s_rho'][::-1]
s_w = variables['s_w'][::-1]

# We transpose (.T) 2D ROMS variables so they have the index ordering expected
# for ISOMIP+ output (ny, nx).  We remove the first and last cell in x and y
# because these are buffer cells in ROMS.
h = variables['h'][1:-1, 1:-1].T
zice = variables['zice'][1:-1, 1:-1].T
Cs_r = variables['Cs_r'][::-1]

Cs_w = variables['Cs_w'][::-1]

mask_rho = variables['mask_rho'][1:-1, 1:-1].T == 1.0
mask_u = variables['mask_v'][1:-1, 1:-1].T == 1.0
mask_v = variables['mask_u'][1:-1, 1:-1].T == 1.0

(nyROMS, nxROMS) = h.shape
nz = len(Cs_r)

hwater = h + zice

# 1D x and y of rho points
x_rho = x0ROMS + dx*(0.5 + numpy.arange(nxROMS))
y_rho = y0ROMS + dy*(0.5 + numpy.arange(nyROMS))

ncOut = Dataset(args.output_file, 'r+')
outVariables = ncOut.variables

for tIndex in range(nTime):
    print "time index {} of {}".format(tIndex, nTime)
    zeta = variables['zeta'][tIndex, 1:-1, 1:-1].T

    # Note: the meaning of u and v are in ROMS are different from
    # u and v in ISOMIP+ because x is lat and y is -lon for the latter
    u = numpy.transpose(variables['v'][tIndex, ::-1, 1:-1, 1:-1], (0, 2, 1))
    v = numpy.transpose(-variables['u'][tIndex, ::-1, 1:-1, 1:-1], (0, 2, 1))

    # -----------------------------------------------------------------------
    #  New formulation: Compute vertical depths (meters, negative) at
    #                   RHO- and W-points, and vertical grid thicknesses.
    #  Various stretching functions are possible.
    #
    #         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
    #
    #                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
    #
    # -----------------------------------------------------------------------

    z_w = numpy.zeros((nz+1, nyROMS, nxROMS))
    z_rho = numpy.zeros((nz+1, nyROMS, nxROMS))
    for zIndex in range(nz+1):
        Zo_w = (hc*s_w[zIndex] + Cs_w[zIndex]*hwater)/(hc + hwater)
        z_w[zIndex, :, :] = mask_rho*(zeta + zice + (zeta + hwater)*Zo_w)
    for zIndex in range(nz):
        Zo_rho = (hc*s_rho[zIndex] + Cs_r[zIndex]*hwater)/(hc + hwater)
        z_rho[zIndex, :, :] = mask_rho*(zeta + zice + (zeta + hwater)*Zo_rho)

    Hz = z_w[0:-1, :, :] - z_w[1:, :, :]

    H_u = 0.5*(Hz[:, :, 0:-1] + Hz[:, :, 1:])
    H_v = 0.5*(Hz[:, 0:-1, :] + Hz[:, 1:, :])

    uMask = numpy.zeros((nz, nyROMS, nxROMS-1), bool)
    vMask = numpy.zeros((nz, nyROMS-1, nxROMS), bool)

    for zIndex in range(nz):
        uMask[zIndex, :, :] = mask_u
        vMask[zIndex, :, :] = mask_v

    # locations of corners in the x-z plane (used for plotting with pcolor and
    # for computing areas of faces involved in momentum fluxes)
    zInterface_u = numpy.zeros((nz+1, nyROMS, nxROMS-1))
    weight = numpy.zeros((nz+1, nyROMS, nxROMS-1))
    for zIndex in range(nz+1):
        zWeight = mask_rho*z_w[zIndex, :, :]
        zInterface_u[zIndex, :, :] = zWeight[:, 0:-1] + zWeight[:, 1:]
        weight[zIndex, :, :] = 1.0*mask_rho[:, 0:-1] + 1.0*mask_rho[:, 1:]
    mask = weight > 0.
    zInterface_u[mask] /= weight[mask]

    osfOut = numpy.ma.masked_all((nzOut, nxOut))
    result = computeOSF(u, uMask, dy=dy*numpy.ones(u.shape),
                        zInterface=zInterface_u,
                        zInterfaceOut=zInterfaceOut, plot=args.plot,
                        xPlot=x_rho, zPlot=z_w, tIndex=tIndex)
    osfOut[:, 1:-1] = result[nzExtra:, :]

    outVariables['overturningStreamfunction'][tIndex, :, :] = osfOut


ncROMS.close()
ncOut.close()
