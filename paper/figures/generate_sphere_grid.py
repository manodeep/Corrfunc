from __future__ import print_function, division

from Corrfunc.utils import gridlink_sphere
import numpy as np
from numpy import sin, cos
from mayavi import mlab


def lcm(arr):
    def gcd(a, b):
        if a == b: return a
        while b > 0: a, b = b, a % b
        return a

    def _lcm(a, b):
        return abs((a // gcd(a, b)) * b)

    return reduce(lambda x, y: _lcm(x, y), arr)    


def spherical_meshgrid(sphere, num_ra, max_dec_cells=200):
    import numpy as np

    ndec = num_ra.size
    grid_ndec_refine = max_dec_cells // ndec + 1
    dec_bands = []
    off = 0
    for idec in xrange(ndec):
        d = sphere['dec_limit'][off]
        dec_bands.extend([d[0]])
        off += num_ra[idec]

    # Take the declination in the last cell; tuple of (d[0], d[1])
    # and then take the upper limit -> d[1] 
    dec_bands.extend([sphere['dec_limit'][-1][1]])

    # Now find the RA max/min
    r = []
    for a in sphere['ra_limit'][0:num_ra[0]]:
        r.extend(a)
        
    ra_limits = [min(r), max(r)]
    dec_cells = np.linspace(dec_bands[0], dec_bands[-1],
                            num=ndec*grid_ndec_refine)
    ra_cells = np.linspace(ra_limits[0], ra_limits[1],
                           num=ndec*grid_ndec_refine)
    dec, ra = np.meshgrid(dec_cells, ra_cells, indexing='ij')
    scalars = np.zeros_like(dec)
    start = 0
    for idec in xrange(ndec):
        scalars[idec*grid_ndec_refine:(idec + 1)*grid_ndec_refine, :] = idec
    
    return dec, ra, dec_bands, scalars


def main(thetamax=None, link_in_ra=True,
         dec_refine_factor=1, ra_refine_factor=1):

    from math import radians
    
    if not thetamax:
        thetamax = 30
    ra_limits = [0.0, 360.0]
    dec_limits = [-90.0, 90.0]
    sphere, num_ra = gridlink_sphere(thetamax,
                                     link_in_ra=link_in_ra,
                                     dec_refine_factor=dec_refine_factor,
                                     ra_refine_factor=ra_refine_factor,
                                     ra_limits=ra_limits,
                                     dec_limits=dec_limits,
                                     return_num_ra_cells=True)

    ra_limits = [radians(a) for a in ra_limits]
    dec_limits = [radians(a) for a in dec_limits]
    
    ndec = num_ra.size
    ncells = sphere.size
    
    dec, ra, dec_bands, scalars = spherical_meshgrid(sphere,
                                                     num_ra)

    x = cos(dec) * cos(ra)
    y = cos(dec) * sin(ra)
    z = sin(dec)

    volume = mlab.mesh(x, y, z, scalars=scalars,
                       colormap='viridis', opacity=0.95)
    volume.actor.property.specular = 0.45
    volume.actor.property.specular_power = 5
    # Backface culling is necessary for more a beautiful transparent
    # rendering.
    volume.actor.property.backface_culling = True


    # Mark the declination bands
    phi = np.linspace(ra_limits[0], ra_limits[1], 100)
    for angle in dec_bands:
        # Lower limits of the dec-bands, except for the two poles
        if angle == -0.5*np.pi or angle == 0.5*np.pi:
            continue
        
        # Draw the declination band
        xx = np.cos(phi) * np.cos(angle)
        yy = np.sin(phi) * np.cos(angle)
        zz = np.ones_like(phi) * np.sin(angle)

        mlab.plot3d(xx, yy, zz, color=(1, 1, 1),
                    opacity=0.5, tube_radius=None)

    # Now draw the ra cells. 
    from itertools import izip, count
    off = 0
    for idec, dec_low, dec_hi in izip(count(), dec_bands[0:-1], dec_bands[1:]):
        dodgerblue = (0.1167315175, 0.5625, 1.0)
        white = (1.0, 1.0, 1.0)
        gray = (0.5, 0.5, 0.5)
        for ira in xrange(num_ra[idec]):
            phi = sphere['ra_limit'][off]
            color, tube_radius = (dodgerblue, 0.01) if ira == 0 else (white, None)
            theta = np.linspace(dec_low, dec_hi, 100)
            xx = np.cos(theta) * np.cos(phi[0]) 
            yy = np.cos(theta) * np.sin(phi[0])
            zz = np.sin(theta)
            mlab.plot3d(xx, yy, zz, color=color,
                        opacity=0.5, tube_radius=tube_radius)
            
            off += 1

        theta = np.linspace(dec_low, dec_hi, 100)
        xx = np.cos(theta) * np.cos(phi[1]) 
        yy = np.cos(theta) * np.sin(phi[1])
        zz = np.sin(theta)
        mlab.plot3d(xx, yy, zz, color=dodgerblue,
                    opacity=0.5, tube_radius=0.01)
        
        
    mlab.show()
    


if __name__ == '__main__':
    main()

    
    
