#############################################################################################
THIS IS OUTDATED IN FILE NAMES AND NUMBERS BUT IS THE GENERAL IDEA OF HOW TO CLEAN UP THE OSC 
#############################################################################################

file: sne_sample.csv

This table is what's left of the Open Supernova Catalog after filtering through entries. These filters include keeping all entries that have discovery dates, host names, and known PGC numbers. The table serves as the candidate catalog to cross-match with the atlas.

HOW TO READ THE TABLE:

    NAME, DATE, HOST, PGC, RA, DEC, TYPE, mmax, z
    NAME: supernova name
    DATE: supernova discovery date year
    PGC: host galaxy pgc number
    HOST: host galaxy name/alias
    RA: supernova right ascension in degrees
    DEC: supernova declination in degrees
    TYPE: supernova type
    mmax: maximum apparent magnitude (?)
    z: redshift
    
HOW TO BUILD THE TABLE (OLD VERSION):

    1.) Download the Open Supernova Catalog (columns: Name, Disc. Date, Host, mmax, RA, Dec, Type, z) as a csv file
    
    2.) clean_catalog.ipynb
    
       Input:
      -Read in the entire Open Supernova Catalog
      -Create function fixed_catalog that does the following:
          -Average the RAs and DECs and convert them into degrees.
          -Collect supernova names, host galaxy names, discovery dates, RAs, DECs, SNe types, mmaxs, and 
      -Build a table called fixed_cat that contains the information that was just collected in the previous function.
      -Read in catalog_pgcs_beta.csv, which is the output of running the catalog through galbase and finding the corresponding PGC number while also putting -1 next to those without a PGC number. This creates catalog_pgcs.csv.
      -Merge fixed_cat with catalog_pgcs, which is a file that was run through galbase in order to see which supernovae have host galaxies with identifiable PGC numbers. This creates catalog_pgcs_total.csv.
      -Loop over catalog_pgcs_total.csv:
          -Get rid of entries without discovery dates. These entries are believed to be supernova remnants.
          -Get rid of entries without host names.
          -Collect all other information.
      -Create table as cleanest_catalog.csv.
      
      Output:
      -cleanest_catalog.csv
      -columns: NAME, DATE, PGC, HOST, RA, DEC, TYPE, mmax, z
      
UPDATE: cleanest_catalog_3.0.csv
    - pretty much the same but with more recent version of the catalog
    - had to rebuild clean_catalog (now 'rebuild clean_catalog')
    - fixed_catalog.csv is the aftermath of getting rid of entries without dates and converting ras and decs to degrees after averaging
    - fixed_catalog_pgcs.csv is fixed catalog but with the added PGCs column (including the flags)
    - cleanest_catalog_3.0.csv is the fixed catalog but minus the entries with -1 for PGC (no recognizable host name)
    
    
HOW TO BUILD THE TABLE (UPDATED VERSION):

    1.) Download the OSC.
    
    2.) sne_sample.ipynb
    
        -Read in the full OSC
        -Create a function fix_catalog that:
            -Averages the RAs and DECs and converts them to degrees
            -Gets rid of entries without RAs or DECs
            -Gets rid of entries without discovery dates (aka get rid of SNRs)
            -Converts all discovery dates to just years
        -Run this function and build the table using pandas DataFrames. Columns: NAME, DATE, HOST, RA, DEC, TYPE, mmax, z
        -Save the table as clean_osc.csv.
        -Run the catalog against galbase and extract PGC info based on host name. If there is no recognizable host name, flag with a '-1'. NOTE: the Python version of galbase is super slow. If you have someone else redo this, then you can skip this step and just read in the file that they send.
        -Add the PGCs column to clean_osc.csv. Might need to do some manipulation with data types (we want strings). Save this file as clean_osc_plus_pgcs.csv.
        -Read in clean_osc_plus_pgs.csv and get rid of all entries with '-1' in the PGC column.
        -Save this as usable_osc.csv. (Number of entries: 2057).