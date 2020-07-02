import numpy
import matplotlib.pyplot as plt
import os


def computeOSF(u, uMask, dy, zInterface, zInterfaceOut, plot=False, xPlot=None,
               zPlot=None, tIndex=None):
    """
    Compute the overturning streamfunction.

    Arguments:
    u -- velocity in x (north) direction with shape (nz, ny, nx-1)

    uMask -- mask of where u is valid with the same shape as u

    dy -- the size of each cell in the y (east-west) with the same shape as u

    zInterface -- the location of interfaces between vertical layers at the
                  same horizontal locations as u, with shape (nz+1, ny, nx-1)

    zInterfaceOut -- the location of interfaces between vertical layers on the
                     output (z-level) grid, with shape nzOut+1

    Keword arguments:
    plot -- if True, plots a section of the velocity before and after
            interpolation to the output grid as well as the zonally integrated
            meridional transport and the overturning streamfunction

    xPLot -- locations of cell centers in x, used to plot the results, with
             shape nx

    zPLot -- locations of layer interfaces as horizontal cell centers, used to
             plot the results, with shape (nz+1, ny, nx)

    tIndex -- the time index, used to determine a file name for plot images

    Returns:
        osfOut -- a masked array containing the overturning streamfunciton on
                  the ISOMIP+ grid with shape (nzOut, nx)
    """

    nz = u.shape[0]
    ny = u.shape[1]
    nx = u.shape[2]+1

    nzOut = len(zInterfaceOut)-1

    uTransportOut = numpy.zeros((nzOut, ny, nx-1))
    uTransportMask = numpy.zeros((nzOut, ny, nx-1), int)
    # now, we need to conservatively interpolate uTransport to a z-level grid
    for zIndexOut in range(nzOut):
        for zIndexIn in range(nz):
            # find the overlap (if any) between this layer on the input and
            # output grids
            zTop = numpy.minimum(zInterface[zIndexIn, :, :],
                                 zInterfaceOut[zIndexOut])
            zBot = numpy.maximum(zInterface[zIndexIn+1, :, :],
                                 zInterfaceOut[zIndexOut+1])

            mask = numpy.logical_and(zTop > zBot, uMask[zIndexIn, :, :])
            area = mask*(zTop-zBot)*dy[zIndexIn, :, :]
            uTransportOut[zIndexOut, :, :] += area*u[zIndexIn, :, :]
            uTransportMask[zIndexOut, :, :] += numpy.array(mask, int)

    # the meridional transport is the sum of the transport in individual cells
    # across the y axis (axis=1)
    uTransportMerid = numpy.sum(uTransportOut, axis=1)
    # uTransportMask and uTransportMeridMask actually contain the number of
    # grid cells on the input grid that contributed to each point in
    # uTransportOut and uTransportMerid, respectively
    uTransportMeridMask = numpy.sum(uTransportMask, axis=1)

    # the overturning streamfunction is the cumulative vertical sum of the
    # meridional transport
    osf = numpy.zeros((nzOut+1, nx-1))
    osfMask = numpy.zeros(osf.shape, bool)
    osfMask[1:, :] = uTransportMeridMask > 0

    osf[1:, :] = numpy.cumsum(uTransportMerid, axis=0)

    # ISOMIP+ output grid is at cell centers, whereas osf is at u points
    # horizontally and w points vertically, so it needs to be interpolated.

    # first average in z
    osf_u = osf[0:-1, :]*osfMask[0:-1, :] + osf[1:, :]*osfMask[1:, :]
    weight_u = 1.0*osfMask[0:-1, :] + 1.0*osfMask[1:, :]

    # then eaverage in x keeping in mind that ISOMIP+ has one extra cell in
    # each direction than the ROMS output
    osfOut = numpy.zeros((nzOut, nx))
    weight = numpy.zeros((nzOut, nx))
    osfOut[:, 0:-1] += osf_u
    osfOut[:, 1:] += osf_u
    weight[:, 0:-1] += weight_u
    weight[:, 1:] += weight_u

    mask = weight > 0.
    osfOut[mask] /= weight[mask]
    osfOut = numpy.ma.masked_array(osfOut, mask=numpy.logical_not(mask))

    if plot:

        if not os.path.exists('plots'):
            os.makedirs('plots')

        uOut = numpy.ma.masked_all((nzOut, ny, nx-1))
        for zIndexOut in range(nzOut):
            zTop = numpy.minimum(zInterface[0, :, :],
                                 zInterfaceOut[zIndexOut])
            zBot = numpy.maximum(zInterface[-1, :, :],
                                 zInterfaceOut[zIndexOut+1])
            areaOut = dy[0, :, :]*(zTop-zBot)
            mask = areaOut > 0.

            uSlice = numpy.zeros((ny, nx-1))
            uSlice[mask] = uTransportOut[zIndexOut, :, :][mask]/areaOut[mask]
            uOut[zIndexOut, :, :] = uSlice

        yIndex = ny/2

        XIn = xPlot.reshape((1, nx))
        ZIn = zPlot[:, yIndex, :]

        plt.close('all')

        plt.figure(1)
        mask = numpy.logical_not(uMask[:, yIndex, :])
        plt.pcolor(1e-3*XIn, ZIn, numpy.ma.masked_array(u[:, yIndex, :],
                                                        mask=mask))
        plt.colorbar()
        plt.title('u (m/s) on the input grid')
        plt.xlabel('x (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/uIn_{:03d}.png'.format(tIndex))

        (XOut, ZOut) = numpy.meshgrid(xPlot, zInterfaceOut, indexing='xy')

        plt.figure(2)
        mask = uOut[:, yIndex, :] == 0.
        plt.pcolor(1e-3*XOut, ZOut, numpy.ma.masked_array(uOut[:, yIndex, :],
                                                          mask=mask))
        plt.colorbar()
        plt.title('u (m/s) on the output grid')
        plt.xlabel('x (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/uOut_{:03d}.png'.format(tIndex))

        plt.figure(3)
        mask = uTransportMeridMask == 0
        plt.pcolor(1e-3*XOut, ZOut, 1e-6*numpy.ma.masked_array(uTransportMerid,
                                                               mask=mask))
        plt.colorbar()
        plt.title('meridional transport (Sv) on the output grid')
        plt.xlabel('x (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/meridTransOut_{:03d}.png'.format(tIndex))

        plt.figure(4)
        plt.pcolor(1e-3*XOut, ZOut, uTransportMeridMask)
        plt.colorbar()
        plt.title('input grid points contributing to meridional tranport')
        plt.xlabel('x (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/meridTransCount_{:03d}.png'.format(tIndex))

        zOut = 0.5*(zInterfaceOut[0:-1] + zInterfaceOut[1:])
        (XOut, ZOut) = numpy.meshgrid(xPlot, zOut, indexing='xy')

        plt.figure(5)
        mask = numpy.logical_not(osfMask[1:-1, :])
        plt.pcolor(1e-3*XOut, ZOut, 1e-6*numpy.ma.masked_array(osf[1:-1, :],
                                                               mask=mask))
        plt.colorbar()
        plt.title('meridional overturning streamfunction (Sv) on the output '
                  'grid')
        plt.xlabel('x (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/osf_{:03d}.png'.format(tIndex))

        dx = xPlot[1]-xPlot[0]
        x0 = xPlot[0]-0.5*dx
        (XOut, ZOut) = numpy.meshgrid(x0 + dx*(numpy.arange(nx+1)),
                                      zInterfaceOut,
                                      indexing='xy')

        plt.figure(6)
        plt.pcolor(1e-3*XOut, ZOut, 1e-6*osfOut)
        plt.colorbar()
        plt.title('meridional overturning streamfunction (Sv) ISOMIP+ grid')
        plt.xlabel('x (km)')
        plt.ylabel('z (m)')
        plt.axis('tight')
        plt.savefig('plots/osfInterp_{:03d}.png'.format(tIndex))

    return osfOut
