---------------------------------------------------
 The Local Environments of Low-Redshift Supernovae
---------------------------------------------------

A project that locates supernovae in nearby galaxies and characterizes
the dust and stellar properties at each explosion site.

################################################################
mask_interp/ and image_processing.ipynb
################################################################

Building of the convolved images were done outside this project (see the .py files in mask_interp/). This includes masking each WISE and GALEX image, blaking foreground stars and galaxies, interpolating the images to fill in these blanked pixels, and convolving them to 2-kpc resolution.
    
We then build the rgrid maps for each galaxy and reproject the FUV maps onto W4 (image_processing.ipynb)

################################################################
sne_sample.ipynb
################################################################

Running a large catalog against galbase can be time-consuming. We grab the PGCs, hosts, inclinations, t values, and tags of galaxies in both the z0MGS sample and galbase. We make csv files with this information for both samples.

    * galbase_info.csv
    
The supernova sample comes from The Open Supernova Catalog (Guillochon et al. 2017).
We take the OSC and clean it up so that it only includes known, historical supernovae with averaged RAs and Decs.

We then crossmatch this version of the OSC with galbase_info.csv to add PGCs, inclinations, and t values. We drop the tags column and output to a csv file. This is our completely cleaned-up version of the OSC that has recognizable host names.

    * clean_osc.csv [only includes types Ia, II, and Ib/c; no IIn or Ia-Pec]
    * clean_osc_full.csv [includes all SN types]

We cross-match the cleaned-up version of the OSC and the z0MGS atlas (Leroy et al. 2019).

    * xmatch.csv [only includes types Ia, II, and Ib/c; no IIn or Ia-Pec]
    * xmatch_full.csv [includes all SN types]
    
We make sure we are confident about each SN by adding columns to the table about host galaxy redshift, the difference between host galaxy redshift and SN redshift, and a believable vs. doubtful sample. We use 2*R25 to indicate whether or not a supernova is inside a galaxy. We target duplicates this way. If duplicates still arise, we check by eye. We also indicate a "believable" sample (i.e.: has either redshift, a host name, or both) and a "doubtful" sample (has neither redshift nor a host name). We create columns that contain host galaxy recessional velocity, host galaxy redshift, and the difference between the redshift of the galaxies and the redshift of the SNe. We exclude galaxies without a position angle or inclination.

    * potential_sample.csv [only includes types Ia, II, and Ib/c; no IIn or Ia-Pec]
    * potential_sample_full.csv [includes all SN types]
    
We target any remaining duplicates by checking by eye.

    * sne_sample.csv [only includes types Ia, II, and Ib/c; no IIn or Ia-Pec]
    * sne_sample_full.csv [includes all SN types]

################################################################
sample_processing.ipynb
################################################################

We build control tables for every galaxy in sne_sample.csv. For every pixel in each galaxy, we measure each WISE and GALEX flux, calculate Sigma_SFR(FUV+W4) and Sigma_SFR(NUV+W3), and compute the galactocentric radius.

    * galaxy_control_tables/

We then use these control tables to build host galaxy CDFs for each band and SFR tracer.

    * galaxy_control_tables/gal_cdfs_2kpc/
    
We next add more information to the SN table, including galactocentric radius, WISE and GALEX fluxes, SFR values, and RMS values.

    * 2kpc.csv
    
Cross-check that these are indeed in 2-kpc resolution images and mark them with a "1" if so.

    * final_2kpc_sample.csv

Using the z0MGS delivery table, we also add information about each SN's' host galaxy, including log10(M_star/M_sun), log10(SFR), and distance. This is all output to a file.

    * FINAL_SAMPLE.csv
    
################################################################
cronin21_sn_catalog.ipynb
################################################################

Create the final SN table in a machine-readable format to then be uploaded as an online catalog.

Read in FINAL_SAMPLE.csv and drop columns to only keep the most relevant ones (see the SN table in the appendix of Cronin et al. 2021). Re-name the columns to have units. Save this as an enhanced csv file (.ecsv).

    * cronin21_sample.ecsv

################################################################
sfms.ipynb
################################################################

We plot the SNe in our sample in the Sigma_SFR / Sigma_star plane to see how SN locations compare to the resolved star-forming main sequence. This involves adding Sigma_SFR(NUV+W3) to the galaxy control tables. We also output the Sigma_SFR and Sigma_star values for each pixel within r25 of each host galaxy.

    * sfr_star.csv

We plot these host galaxy values and the values for each SN.

    * sfr_vs_star_full.png [all SN types together]
    * sfr_vs_star.png [splitting the plot into 3 subplots: SNe Ia, II, and Ib/c]
    
We next made a KDE plot to show how each SN distributions as a function of Sigma_SFR and as a function of Sigma_SFR/Sigma_star compare with each other. We quantify this comparison by performing a T-test between each SN type.

    * sfr_kde.pdf
    * sfr_star_kde.pdf
    
We then combine SNe II and Ib/c into one category of CC SNe and perform the same analysis above, this time looking at how the distributions of SNe Ia and CC SNe change by galaxy type (i.e., early- vs. late-type galaxies).

    * sfr_type_kde_v2.pdf
    * sfr_star_type_kde_v2.pdf
    

################################################################
monte_carlo.ipynb
################################################################

Run a Monte Carlo simulation that will randomly sample pixels from galaxies in each band and SFR tracer. Let the sample size of each simulated distribution be equal to the number of SNe per type in each tracer. Repeat this 100 times for 100 total simulated distributions. This will show us the uncertainty/spread in the SN CDFs that is due to smaller sample sizes. Output each 100 distribution to a single file per SN type per band.

    * monte_carlo/

################################################################
results.ipynb
################################################################

Build the SN cdfs for each WISE and GALEX band and SFR tracer. Calculate toy models of distributions that are 2x more extended and 2x more concentrated than host galaxy emission. Add the Monte Carlo distributions to the plots to visualize the uncertainty.

    * final_results_IR.pdf
    * final_results_UVSFR.pdf
    
Output the CDFs.

    * 2kpc.csv
    * SFR_2kpc.csv
    
Calculate KS and Anderson-Darling tests using the Monte Carlo distributions and the SN distributions.

    * ks_2kpc_mc.csv [using Monte Carlo p-values; this is a combination of the above files]
    * ad_2kpc_mc.csv [Anderson-Darling version]
    
    
We also build some extra plots that help the methods section of the paper. We first plot some sample W1 images of galaxies and mark their SNe.

    * fig1.png
    
We also plot a sample host galaxy CDF with its SNe marked.

    * gal_cdf.pdf
    
We then plot the toy models of distributions that are 2x more extended and 2x more concentrated than host galaxy emission and show how they change when going from the host galaxy CDF to the SN CDFs.

    * fiducial_flux.pdf
    * fiducial_gal_cdf.pdf
