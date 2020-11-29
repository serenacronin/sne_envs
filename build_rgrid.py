from deproject import deproject
import numpy as np
import pandas as pd
from astropy.wcs import WCS



def build_rgrid(wcs, ra_gal, dec_gal, incl_gal, pa_gal):
    """ Input: WCS of the fits image for a galaxy
        Output: a grid of galactocentric radii for each pixel in the galaxy image
    """
    naxis = wcs._naxis # grab the axes
    
    # make an map 'x' for pixels along the RA axis and a map 'y' for pixels along the DEC axis
    x, y = np.meshgrid(range(naxis[0]), range(naxis[1]))
    
    # convert each x and y map to RA and DEC maps in degrees
    ra_deg, dec_deg = wcs.wcs_pix2world(x,y, 0)
    
    rgal_map, phigal_map = deproject(center_coord=(ra_gal, dec_gal), incl=incl_gal, pa=pa_gal, wcs=wcs, naxis=naxis)
    
    return(rgal_map)
