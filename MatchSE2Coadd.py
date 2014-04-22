#!/usr/bin/env python

import sys
import os
import subprocess
import pyfits
import numpy as np

from esutil import htm
import warnings


def get_from_list( listfile, prefix ):
    list = open(listfile).read().strip().split()
    for i in range(len(list)):
        list[i] = '%s%s' %(prefix,list[i])
    return list

def not_unique(arr):
    sorted = np.sort(arr)
    diff = np.diff(sorted)
    diff = np.r_[1,diff]

    not_unique = np.unique( sorted[ diff==0 ] )
    #not_unique = sorted[ diff==0 ]
    return not_unique

def to_delete(from_m, to_m, rad, not_unique):
    nq = (to_m==not_unique)
    indecies = np.where( nq==True )[0]
    dups = from_m[nq]
    seps = rad[nq]
    ind = np.argmin(seps)
    dup_index = indecies[ind]
    to_delete = indecies[ indecies!=dup_index ]
    return to_delete

def scan_dup(to_be_deleted, from_m, to_m, rad, not_unique):
    for i in range(len(not_unique)):
        dels = to_delete(from_m, to_m, rad, not_unique[i])
        to_be_deleted = np.append( to_be_deleted, dels )
    return to_be_deleted

def one_to_one_closest(m1, m2, radius):
    to_be_deleted = np.array( [] )

    not_uni = not_unique(m1)
    to_be_deleted = scan_dup( to_be_deleted, m2, m1, radius, not_uni )

    not_uni = not_unique(m2)
    to_be_deleted = scan_dup( to_be_deleted, m1, m2, radius, not_uni )

    m1 = np.delete(m1,to_be_deleted)
    m2 = np.delete(m2,to_be_deleted)
    radius = np.delete(radius,to_be_deleted)

    return [m1,m2,radius]




def cat_catalogs( exp_list, se_dir, allchips_dir, se_suffix ):
    cats = []
    for exp in exp_list:
        exp = exp.rstrip('/')

        if allchips_dir == None:
            exp_objs_file = os.path.join( se_dir, exp, '%s_allchips.fits'%(exp) )
        else:
            exp_objs_file = os.path.join( allchips_dir, '%s_allchips.fits'%(exp) )
        cats.append(exp_objs_file)

        if os.path.exists(exp_objs_file):
            subprocess.call( ['rm', exp_objs_file] )

        files = []
        for i in range(1, 63):
            if str(i) in skip_ccd:
                continue

            if i < 10:
                index = '0%i' %i
            else:
                index = '%i' %i
            file = '%s_%s%s' %(exp,index,se_suffix)
            file = os.path.join(se_dir, exp, file)
            files.append(file)

        file = files[0]
        hdus = pyfits.open(file)
        tbhdu = hdus[2]
        names = tbhdu.data.columns.names
        formats = tbhdu.data.columns.formats

        exp_objs = tbhdu.data
        for file in files[1:]:
            hdus = pyfits.open(file)
            exp_objs = np.hstack( (exp_objs, hdus[2].data) )
        
        cols = []
        for i in range(len(names)):
            c = pyfits.Column( name=names[i], format=formats[i], array=exp_objs[names[i]] )
            cols.append(c)
        c = pyfits.Column( name='INDEX', format='J', array=np.arange(len(exp_objs)) )
        cols.append(c)

        columns = pyfits.ColDefs(cols)
        tbhdu = pyfits.new_table(columns)
        phdu = pyfits.PrimaryHDU()
        hdus = pyfits.HDUList( [phdu,tbhdu] )
        if os.path.exists(exp_objs_file):
            subprocess.call( ['rm', exp_objs_file] )
        hdus.writeto(exp_objs_file)
    return cats


if __name__=='__main__':
    warnings.filterwarnings('ignore')

    explist = sys.argv[1]
    coadd_cat = sys.argv[2]
    se_dir = sys.argv[3]
    outfile = sys.argv[4]
    atleast = int(sys.argv[5])
    allchips_dir = sys.argv[6]
    se_suffix = sys.argv[7]
    prefix = sys.argv[8]
    skip_ccd = sys.argv[9]
    skip_ccd = skip_ccd.split(',')

    if allchips_dir == 'None':
        allchips_dir = None

    coadd_hdus = pyfits.open( coadd_cat )
    coadd_hdu = coadd_hdus[1]
    coadd_data = coadd_hdu.data
    coadd_ra = coadd_data['ALPHAWIN_J2000']
    coadd_dec = coadd_data['DELTAWIN_J2000']

    exp_list = get_from_list(explist, prefix)
    combined = cat_catalogs( exp_list, se_dir, allchips_dir, se_suffix )
    cl = len(combined)

    h = htm.HTM()
    maxrad = 1.0/3600.0

    se_ras = np.array( [-999.]*len(coadd_data)*cl).reshape(len(coadd_data),cl)
    se_decs = np.array( [-999.]*len(coadd_data)*cl).reshape(len(coadd_data),cl)
    se_indecies = np.array( [-1.]*len(coadd_data)*cl).reshape(len(coadd_data),cl)

    for i in range(len(combined)):
        
        se_file = combined[i]
        se_data = pyfits.open( se_file )[1].data
        se_ra = se_data['ALPHAWIN_J2000']
        se_dec = se_data['DELTAWIN_J2000']

        m1,m2,radius = h.match(se_ra,se_dec,coadd_ra,coadd_dec,maxrad,maxmatch=1)
        m1,m2,radius = one_to_one_closest(m1,m2,radius) 
        
        for j in range(len(m1)):
            se_ras[ m2[j] ][i] = se_ra[ m1[j] ] 
            se_decs[ m2[j] ][i] = se_dec[ m1[j] ]
            se_indecies[ m2[j] ][i] = m1[j]
   
    count = 0
    for i in range(len(se_indecies)):
        line_cut = ( se_indecies[i]!=-1 )
        ok = np.sum(line_cut)
    
    cols = []
    coadd_names = coadd_hdu.columns.names
    coadd_formats = coadd_hdu.columns.formats
    for name, format in zip(coadd_names,coadd_formats):
        col = pyfits.Column( name=name, format=format, array=coadd_data[name] )
        cols.append(col)
        

    se_names = ['SE_INDEX','SE_RA','SE_DEC']
    se_formats = ['%iJ'%cl, '%iD'%cl, '%iD'%cl, '%iD'%cl]
    se_arrays = [se_indecies, se_ras, se_decs]
    for name,format,array in zip(se_names,se_formats,se_arrays):
        col = pyfits.Column( name=name, format=format, array=array )
        cols.append(col)

    columns = pyfits.ColDefs(cols)
    tbhdu = pyfits.new_table(columns)
    phdu = pyfits.PrimaryHDU()
    hdus = pyfits.HDUList( [phdu,tbhdu] )
    if os.path.exists(outfile):
        subprocess.call( ['rm', outfile] )
    hdus.writeto(outfile)
    
