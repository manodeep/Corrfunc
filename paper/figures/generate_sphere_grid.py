from __future__ import print_function, division

from Corrfunc.utils import gridlink_sphere
import numpy as np
from numpy import sin, cos
from mayavi import mlab


# def gcd(*numbers):
#     """Return the greatest common divisor of the given integers"""
#     from functools import reduce
#     from fractions import gcd
#     return reduce(gcd, numbers)

# # Least common multiple is not in standard libraries? It's in gmpy, but this is simple enough:
# def lcm(*numbers):
#     """Return lowest common multiple."""
#     from functools import reduce

#     def lcm(a, b):
#         return (a * b) // gcd(a, b)
#     return reduce(lcm, numbers, 1)




# def lcm(arr):
#     from functools import reduce    # need this line if you're using Python3.x
#     def _lcm(a, b):
#         if a > b:
#             greater = a
#         else:
#             greater = b

#         while True:
#             if greater % a == 0 and greater % b == 0:
#                 lcm = greater
#                 break
#             greater += 1
        
#         return lcm

#     return reduce(lambda x, y: _lcm(x, y), arr)

def lcm(arr):
    def gcd(a, b):
        if a == b: return a
        while b > 0: a, b = b, a % b
        return a

    def _lcm(a, b):
        return abs((a // gcd(a, b)) * b)

    return reduce(lambda x, y: _lcm(x, y), arr)    


def generate_meshgrid_from_sphere(sphere, num_ra_cells,
                                  min_ra_cells_per_band=200):
    
    ndec = num_ra_cells.size
    ncells = sphere.size

    # Find the lowest common multiple of the
    # number of RA cells in any dec band
    ra_cells = lcm(num_ra_cells)

    # But make sure there are at least min_ra_cells_per_band
    fac = long(min_ra_cells_per_band // ra_cells) + 1
    ra_cells *= fac

    assert ra_cells > min_ra_cells_per_band
    
    ## Now create a mesh grid that goes from the min to max
    ## ra and dec
    dec = np.zeros((ndec + 1, ra_cells))
    ra = np.zeros((ndec + 1, ra_cells))
    scalars = np.zeros((ndec + 1, ra_cells))

    dec_start = 0
    for idec in xrange(ndec):
        dec[idec, :] = sphere['dec_limit'][dec_start][0]
        scalars[idec, :] = idec
        refine = ra_cells // num_ra_cells[idec]
        assert refine > 1
        for ira in xrange(num_ra_cells[idec]):
            this_ra = sphere['ra_limit'][dec_start + ira]
            last_cell = True if ira == num_ra_cells[idec] - 1 else False
            
            if last_cell:
                this_refine = refine - 1
            else:
                this_refine = refine

            ra_binsize = (this_ra[1] - this_ra[0])/this_refine
            refined_ra = np.array([this_ra[0] + d*ra_binsize for d in xrange(this_refine)])
            start = ira*refine
            dest_sel = np.s_[start:start + this_refine]
            ra[idec, dest_sel] = refined_ra
            if last_cell:
                ra[idec, -1] = this_ra[1]
            
        dec_start += num_ra_cells[idec]

    dec[ndec, :] = sphere['dec_limit'][dec_start - num_ra_cells[ndec-1]][1]
    scalars[ndec, :] = ndec
    ra[ndec, :] = ra[ndec-1, :]
        
    return dec, ra, scalars


def main(thetamax=30, link_in_ra=True,
         dec_refine_factor=1, ra_refine_factor=1):

    sphere, num_ra = gridlink_sphere(thetamax,
                                     link_in_ra=link_in_ra,
                                     dec_refine_factor=dec_refine_factor,
                                     ra_refine_factor=ra_refine_factor,
                                     return_num_ra_cells=True)

    ndec = num_ra.size
    ncells = sphere.size
    
    dec, ra, scalars = generate_meshgrid_from_sphere(sphere,
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
    phi = np.linspace(0, 2 * np.pi, 100)

    # Draw the equator
    xx = np.cos(phi)
    yy = np.sin(phi)
    zz = np.zeros_like(phi)
    mlab.plot3d(xx, yy, zz, color=(1, 1, 1),
                opacity=1.0, tube_radius=0.01)

    # Lower limits of the dec-bands, except for the two poles
    dec_bands = [d[0] for d in sphere['dec_limit'] if d[0] != -0.5*np.pi  or d[1] != 0.5 * np.pi]
    for angle in dec_bands:

        # Draw the declination band
        xx = np.cos(phi) * np.cos(angle)
        yy = np.sin(phi) * np.cos(angle)
        zz = np.ones_like(phi) * np.sin(angle)

        mlab.plot3d(xx, yy, zz, color=(1, 1, 1),
                    opacity=0.2, tube_radius=None)

    # Now draw the ra cells. Need to extend to the poles

    # Note, includes last element unlike when we
    # generated the declination band-lines
    dec_bands = []
    off = 0
    for idec in xrange(ndec):
        d = sphere['dec_limit'][off]
        dec_bands.extend([d[0]])
        off += num_ra[idec]

    # Take the declination in the last cell; tuple of (d[0], d[1])
    # -> then take d[1]
    dec_bands.extend([sphere['dec_limit'][-1][1]]) 
    from itertools import izip, count

    off = 0
    for idec, dec_low, dec_hi in izip(count(), dec_bands[0:-1], dec_bands[1:]):
        for ira in xrange(num_ra[idec]):
            phi = sphere['ra_limit'][off][0]
            theta = np.linspace(dec_low, dec_hi, 100)
            xx = np.cos(theta) * np.cos(phi) 
            yy = np.cos(theta) * np.sin(phi)
            zz = np.sin(theta)
            mlab.plot3d(xx, yy, zz, color=(1, 1, 1),
                        opacity=0.5, tube_radius=None)
            
            off += 1
    mlab.show()
    


if __name__ == '__main__':
    main()

# def smooth_ra_dec(this_ra, this_dec,
#                   RA_subdivide_factor=10):

#     if RA_subdivide_factor == 1:
#         return this_ra, this_dec
    
#     smoothed_dec = []
#     smoothed_ra = []
#     for jj in xrange(len(this_ra) - 1):
#         next_point = jj + 1
#         if this_ra[jj] == this_ra[next_point]:
#             smoothed_ra.extend([this_ra[jj]])
#             smoothed_dec.extend([this_dec[jj]])
#             continue

#         ## So the connectivity, as specified by this_ra
#         ## is spanning a range of RA. Let's break up into
#         ## smaller RA cells
        
#         # But first make sure that the dec's are the same
#         assert this_dec[jj] == this_dec[next_point]
        
#         new_ra_binsize = (this_ra[next_point] - this_ra[jj])/ RA_subdivide_factor
#         new_ra_points = [this_ra[jj] + d * new_ra_binsize for d in xrange(RA_subdivide_factor)]
#         assert len(new_ra_points) == RA_subdivide_factor
#         new_dec_points = [this_dec[jj]]*len(new_ra_points)
#         smoothed_dec.extend(new_dec_points)
#         smoothed_ra.extend(new_ra_points)
        
#     smoothed_dec.extend([this_dec[-1]])
#     smoothed_ra.extend([this_ra[-1]])

#     return smoothed_ra, smoothed_dec

    

# def mlab_plot3d_impl(sphere, num_ra_cells, RA_subdivide_factor=10):

#     ndec = num_ra.size
#     ncells = sphere.size
#     assert num_ra.sum() == ncells

#     cellid = np.arange(ncells)

#     grid_dec = []
#     grid_ra = []

#     # vertices are stored in the connectivity order
#     # final shape is going to be 
#     vertices = []
#     vert_cellids = []
#     last_ra_cell = np.zeros(ncells, dtype=np.int64)
#     dec_band_for_cell = np.zeros(ncells, dtype=np.int64)
#     start = 0
#     for idec in xrange(ndec):
#         dest_sel = np.s_[start:start+num_ra[idec]]
#         last_ra_cell[dest_sel] = start + num_ra[idec] - 1
#         dec_band_for_cell[dest_sel] = idec
#         start += num_ra[idec]

#     first_ra_cell_in_dec_band = np.empty_like(num_ra)
#     first_ra_cell_in_dec_band[0] = 0
#     first_ra_cell_in_dec_band[1:] = num_ra.cumsum()[0:-1]

#     dec_band = -1 # hack to make the increment work for all 
#     for ii, d, r in zip(cellid, sphere['dec_limit'], sphere['ra_limit']):

#         # is this the first RA cell in the declination band?
#         dec_band = dec_band_for_cell[ii]
#         scalar = first_ra_cell_in_dec_band[dec_band]

#         # Specify the connectivity for this cell
#         # clockwise
#         this_dec = [d[0], d[0], d[1], d[1], d[0]]
#         this_ra  = [r[0], r[1], r[1], r[0], r[0]]

#         next_point = ii + 1
#         if ii != last_ra_cell[ii]:
#             ## If not the last RA cell, then bring back to the
#             ## starting point for the next RA cell
#             this_dec.extend([d[0]])
#             this_ra.extend([r[1]])
#         else:
#             ## Last RA cell -> bring to next dec band, first RA cell
#             if ii == ncells - 1:
#                 continue

#             this_dec.extend([d[1]])
#             this_ra.extend([r[1]])


#         # Now increase the number of cells by smoothing the RA cell by 
#         # RA_subdivide_factor. This is simply for a better plot
#         smoothed_ra, smoothed_dec = smooth_ra_dec(this_ra, this_dec,
#                                                   RA_subdivide_factor)
#         assert len(smoothed_dec) == len(smoothed_ra)
        
#         grid_dec.extend(smoothed_dec)
#         grid_ra.extend(smoothed_ra)
#         vert_cellids.extend([scalar]*len(smoothed_dec))
#         assert len(smoothed_dec) == len(smoothed_ra)
#         assert len(smoothed_ra) > RA_subdivide_factor
#         assert len(vert_cellids) == len(grid_dec)

#         this_vertices = [(dd, rr) for dd,rr in zip(smoothed_dec, smoothed_ra)]
#         vertices.extend(this_vertices)
#         print("{0:3d}: RA = {1} verts = {2}".format(dec_band, r, this_vertices))
        
#         assert len(grid_dec) == len(vertices)


#     assert len(grid_dec) == len(vertices)
#     assert len(grid_dec) == len(vert_cellids)
#     assert len(grid_dec) == len(grid_ra)
    
#     grid_dec = np.array(grid_dec)
#     grid_ra = np.array(grid_ra)
#     vertices = np.array(vertices)
#     vert_cellids = np.array(vert_cellids)
    
#     x = cos(grid_dec) * cos(grid_ra)
#     y = cos(grid_dec) * sin(grid_ra)
#     z = sin(grid_dec)

#     mlab.plot3d(x, y, z, vert_cellids, colormap='viridis')
#     mlab.show()
    
    
