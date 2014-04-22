#!/usr/bin/env python

import sys
import os
import subprocess


if __name__=='__main__':

    tile = 'DES0102-4914'
    band = 'i'
    compto_band = 'r'
    topdir_coadd = 'desar2.cosmology.illinois.edu/DESFiles/desardata/OPS/coadd/20130805000011_DES0102-4914/coadd/' #where coadd cats live
    outdir = '%s_prererun' %tile #where to write output

    # plotting options
    class_star = 0.98
   

    # You might need to change these
    coadd_suffix = '_cat.fits'
    coadd_cat = os.path.join(topdir_coadd, '%s_%s%s'%(tile,band,coadd_suffix))
    compto_cat = os.path.join(topdir_coadd, '%s_%s%s'%(tile,compto_band,coadd_suffix))



    ##############################################################################################
    ### shouldn't really have to edit this unless you move the script files somewhere else

    exc = os.path.join('./PlotTwoBands.py')

    if not os.path.lexists(outdir):
        os.makedirs(outdir)

    title = '%s_to_%s'%(band,compto_band)
    subprocess.call( [exc, coadd_cat, compto_cat, title, tile, outdir, str(class_star)] )

    
