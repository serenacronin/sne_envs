# INTERPOLATE GALAXY IMAGE AFTER APPLYING STAR AND GALAXY MASKS
# BINS GALAXY BY RADIUS OUT TO THREE TIMES R25 USING A BIN WIDTH EQUAL TO THE ANGULAR RESOLUTION
# INTERPOLATES WITH NP.INTERP USING THE MEDIAN VALUE OF EACH BIN AS A GUIDE

from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.stats import median_absolute_deviation as mad
import astropy.units as u
#from make_rgrid import gen_rgrid
from high_incl_interp import high_incl_interp
from reproject import reproject_interp


# IMAGE AND TABLE PATHS
im_path = '/data/kant/0/leroy.42/allsky/delivery/'            # gal images
mask_path = '/data/kant/0/leroy.42/allsky/masks/'             # gal masks
d_path = '/data/tycho/0/behrens.52/data/gauss'                # data table that contains all PGCs
#interp_path = '/data/tycho/0/behrens.52/masked_interp/'       # interpolated images
interp_path = '/data/kant/0/leroy.42/allsky/postprocess/interpol/'
rgrid_path = '/data/kant/0/leroy.42/allsky/postprocess/new_rgrid/'

bands = ['w1','w2','w3','w4','fuv','nuv']
res_vals = ['15']

for res in res_vals:
    data = fits.getdata(d_path + res + '_3.5kpc.fits')
    pgcs = data['PGC_NUM']
    
    print('\nResolution: ' + res + ' arcsec')

    for band in bands:
        print('\n\tBand: ' + band)
        
        if res == '7p5' and band == 'w4':            # no w4 w/ 7.5" arcsec res
            continue

        for i in range(len(pgcs)):

            if i == int(len(pgcs)/4):
                print('\t\t25% of the way through ' + band)
            elif i == int(len(pgcs)/2):
                print('\t\t50% of the way through ' + band)
            elif i == 3 * int(len(pgcs)/4):
                print('\t\t75% of the way through ' + band)
            
            pgc = pgcs[i]
            
            if pgc == 2557 or pgc == 5818:
                continue

            if os.path.isfile('/data/kant/0/leroy.42/allsky/postprocess/interpol/PGC%s_%s_gauss15_interp.fits' % (pgc,band)) is True:
                continue
            
	    # Defining path and images name
            f = im_path + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '.fits'                         # gal image
            f_star_mask = mask_path + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '_mask_stars.fits'  # star mask
            f_gal_mask = mask_path + 'PGC' + str(pgc) + '_gauss' + res + '_galaxies.fits'                  # galaxy mask
            #f_rgrid = mask_path + 'PGC' + str(pgc) + '_gauss' + res +'_rgrid.fits'                         # rgrid (degrees away from center)
            f_rgrid = rgrid_path  + 'PGC' + str(pgc) + '_gauss15_rgrid.fits'
      
            if not os.path.isfile(f):
                continue
            
            # APPLY STAR AND GAL MASKS AND CHANGE 1'S TO NAN VALS
            im = fits.getdata(f)
            #star_mask = fits.getdata(f_star_mask)
            #gal_mask = fits.getdata(f_gal_mask)
            hdul_star = fits.open(f_star_mask)
            hdul_gal = fits.open(f_gal_mask)
            hdul_im = fits.open(f)
             
            # use this and the above hdul lines for 15 arcsec images
            # reprojecting the stars and galaxies masks to the header of the 15 arcsec images
            star_mask, footprint_stars = reproject_interp(hdul_star, hdul_im[0].header,order='nearest-neighbor')
            gal_mask, footprint_gals,   = reproject_interp(hdul_gal, hdul_im[0].header,order='nearest-neighbor')
                       

            if res == '15' and np.isnan(data['POSANG'][i]):
                continue
            
            rad_deg = data['RAD_DEG'][i]          # R25 in degrees
            rgrid = fits.getdata(f_rgrid)  
            
            if im.shape != rgrid.shape:
                continue
            
            inds = np.where((star_mask == 1) | (gal_mask == 1))
            if len(inds[0]) == 0:      # if nothing is masked within 3 gal radii of gal center, skip galaxy
                continue
           
            im[inds] = np.nan
            masked = np.copy(im)
        
            if res == '7p5':
                step = 7.5*u.arcsec.to(u.degree)                # bin width, degrees
            else:
                step = 15*u.arcsec.to(u.degree)                 # bin width, degrees
          
            # some rgrid maps are full of 'nan's because of missing data (eg. pos_ang)
            try:
                bins = np.arange(rgrid.min(),3*rad_deg+step,step)   # bin galaxy by radius out to 3*R25
            except:
                continue
        
            if bins.size == 0.:
                continue

            med_vals = []

            # GET MEDIAN VALUE OF EACH BIN
            for k in range(len(bins)-1):
                bin_min = bins[k]                                                 # left edge of bin
                bin_max = bins[k+1]                                               # right edge of bin
                bin_inds = np.where((rgrid >= bin_min) & (rgrid <= bin_max))      # indices of pixels in this radius bin
                bin_med = np.nanmedian(im[bin_inds])                              # median value of bin
                
                while np.isnan(bin_med) and (bin_max <= 3*rad_deg):                # if all values in bin are nan vals
                    bin_max += step/2
                    bin_inds = np.where((rgrid >= bin_min) & (rgrid <= bin_max))  # indices of pixels in this radius bin
                    bin_med = np.nanmedian(im[bin_inds])        
                
                med_vals.append(bin_med)

            # INTERPOLATE
            bins_mid = 0.5*(bins[:-1]+bins[1:])
            im[inds] = np.interp(rgrid[inds],bins_mid,med_vals)
           
            # SAVE INTERPOLATED IMAGE
            hdr = fits.getheader(f)
            
            new_f = interp_path + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '_interp.fits'
            #new_f = '/home/cronin.104/Desktop/mask_interp/' + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '_interp.fits'
            fits.writeto(new_f,im,header=hdr,overwrite=True) 

            # PLOT MASKED IMAGE (OR OLD INTERPOLATED IMAGE) NEXT TO NEW INTERPOLATED IMAGE
            #old_im = fits.getdata(interp_path + 'PGC' + str(pgc) + '_' + band + '_gauss' + res + '_interp.fits')
            #incl = data[i]['INCL']
            
            #vmin = np.nanmedian(im) - 5 * mad(im,ignore_nan=True)
            #vmax = np.nanmedian(im) + 20 * mad(im,ignore_nan=True)
            
            #plt.figure(figsize=(10,6))
            #plt.subplot(121)
            #plt.imshow(masked,vmin=vmin,vmax=vmax,origin='lower')
            #plt.imshow(old_im,vmin=vmin,vmax=vmax,origin='lower')
            #plt.subplot(122)
            #plt.imshow(im,vmin=vmin,vmax=vmax,origin='lower')
           # plt.suptitle('PGC ' + str(pgc) + ', Incl = %.2f' % incl,fontsize=20)
            #plt.tight_layout()
            #plt.show()

                    
        
                
   
