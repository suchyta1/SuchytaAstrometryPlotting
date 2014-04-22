#!/usr/bin/env python

import subprocess
import os
import sys


if __name__=='__main__':

    tile = 'DES0102-4914'
    band = 'i'
    topdir_coadd = 'desar2.cosmology.illinois.edu/DESFiles/desardata/OPS/coadd/20130805000011_DES0102-4914/coadd/' #where coadd cats live
    topdir_SE = 'desar2.cosmology.illinois.edu/DESFiles/desardata/OPS/red/20130717213138_20121124/red/' #where SE exp subdirs live
    exposure_list_SE = 'example_i.fits'
    prefix_list_line = 'DECam_00' # This is appended to the beginning of each line in your exposure list files to make file names
    outdir = '%s_prererun' %tile #where to write output

    # plotting options
    match_at_least = 5 # Only objects that show up in this many SE catalogs show up in the plot
    class_star = 0.9 # Only objecst with CLASS_STAR larger than this show up in the plot


    # You might need to change these
    skip_ccd = '61' # Ignore anything to do with chip 61
    SE_suffix = '_cat.fits'
    coadd_suffix = '_cat.fits'
    coadd_cat = os.path.join(topdir_coadd, '%s_%s%s'%(tile,band,coadd_suffix))
    allchips_dir = None # Leaving this as None writes a file into same directory as SE catalogs which concatenates all SE per chip cats into one FP cat (_allchips appended to name). Changing it to something else will write one of these files for each exposure to the given directory.



    ##############################################################################################
    ### shouldn't really have to edit this unless you move the script files somewhere else

    match_SE2coadd_script = './MatchSE2Coadd.py'
    plot_script = './PlotSameBand.py'
 
   
    if not os.path.lexists(outdir):
        os.makedirs(outdir)
    cm = '%s_%s_SE2coadd.fits' %(tile,band)
    cm_file = os.path.join(outdir,cm )

    subprocess.call( [match_SE2coadd_script, exposure_list_SE, coadd_cat, topdir_SE, cm_file, str(match_at_least), str(allchips_dir), SE_suffix, prefix_list_line, str(skip_ccd)] )
    subprocess.call( [plot_script, cm_file, str(match_at_least), tile, band, outdir, str(class_star)] )

