import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.stats import median_absolute_deviation as mad
import matplotlib.pyplot as plt

d_path = '/data/tycho/0/behrens.52/data/gauss' 
im_path = '/data/kant/0/leroy.42/allsky/delivery/'            # gal images
mask_path = '/data/kant/0/leroy.42/allsky/masks/'  

def high_incl_interp(f,index,res,im,rgrid):

    #### REMOVE BELOW BEFORE ADDING FUNCTION TO MASK_INTERPY.PY #### 

    f = im_path + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '.fits'                         # gal image
    f_star_mask = mask_path + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '_mask_stars.fits'  # star mask
    f_gal_mask = mask_path + 'PGC' + str(pgc) + '_gauss' + res + '_galaxies.fits'                  # galaxy mask

    im = fits.getdata(f)
    star_mask = fits.getdata(f_star_mask)
    gal_mask = fits.getdata(f_gal_mask)

    inds = np.where((star_mask == 1) | (gal_mask == 1))

    im[inds] = np.nan
    masked = np.copy(im)

    #### REMOVE ABOVE BEFORE ADDING FUNCTION TO MASK_INTERPY.PY #### 

    data = fits.getdata(d_path + res + '_3.5kpc.fits')
    hdr = fits.getheader(f)
    # refx, refy are coords of point that needs to be interpolated
    # define position angle as pa
    pa = np.radians(data['POSANG'][index])
    theta = pa - np.pi/2   # pa drawn from pos x axis
    rad_deg = data['RAD_DEG'][index]

    # MAKE ARRAYS OF DISTANCES TO MAJOR AND MINOR AXES OF GAL
    # dist = abs(ax + by + c) / sqrt(a**2 + b**2)
    ctr = int(hdr['CRPIX1'])       # center pixel
    maj_slope = np.tan(theta)      # major axis slope, a
    b = -1                      
    c_maj = -maj_slope*ctr + ctr   
    
    min_slope = np.tan(pa)
    c_min = -min_slope*ctr + ctr

    y,x = np.indices(im.shape)
    maj_dist = (abs(maj_slope*x + b*y + c_maj) / np.sqrt(maj_slope**2 + b**2)) * hdr['CD2_2']    # dist to major axis
    min_dist = (abs(min_slope*x + b*y + c_min) / np.sqrt(min_slope**2 + b**2)) * hdr['CD2_2']    # dist to minor axis

    beam = hdr['BMAJ']             # beam in degrees
    strips = np.arange(0,10*beam + 0.5*beam,0.5*beam)

    for i in range(len(strips)):
        
        inds = np.where(((maj_dist >= i*0.5*beam) & (maj_dist <= (i+1)*0.5*beam)) & (min_dist <= 2*rad_deg))                                           # indices of pixels in strip
        nan_inds = np.where(np.isnan(im[inds]))[0]                   
        strip_nanx = inds[1][nan_inds]                               # x values of nan pixels in strip
        strip_nany = inds[0][nan_inds]                               # y values of nan pixels in strip

        bins = np.arange(0,min_dist[inds].max() + 0.5*beam,0.5*beam) # bins from 0 to min_dist = 2*r25

        med_vals = []

        # GET MEDIAN VALUE OF EACH BIN IN STRIP
        for k in range(len(bins)-1):
            bin_min = bins[k]                                                 # left edge of bin
            bin_max = bins[k+1]                                               # right edge of bin
            bin_inds = np.where((min_dist[inds] >= bin_min) & (min_dist[inds] <= bin_max))[0]      # indices of pixels in this radius bin
            binx = inds[1][bin_inds]
            biny = inds[0][bin_inds]
            bin_med = np.nanmedian(im[biny,binx])                              # median value of bin        
            med_vals.append(bin_med)
    
        bins_mid = 0.5*(bins[:-1]+bins[1:])
        im[strip_nany,strip_nanx] = np.interp(min_dist[strip_nany,strip_nanx],bins_mid,med_vals)
    
    incl = data['INCL'][index]

    vmin = np.nanmedian(im) - 5 * mad(im,ignore_nan=True)
    vmax = np.nanmedian(im) + 30 * mad(im,ignore_nan=True)

    plt.figure(figsize=(10,6))
    plt.subplot(121)
    plt.imshow(masked,vmin=vmin,vmax=vmax,origin='lower')
    #plt.imshow(old_im,vmin=vmin,vmax=vmax,origin='lower')
    plt.subplot(122)
    plt.imshow(im,vmin=vmin,vmax=vmax,origin='lower')
    plt.suptitle('PGC ' + str(pgc) + ', Incl = %.2f' % incl,fontsize=20)
    plt.tight_layout()
    plt.show()


