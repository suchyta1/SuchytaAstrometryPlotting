#!/usr/bin/env python

import sys
import os
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt

from esutil import htm
from pylab import *
from util import *


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



if __name__=='__main__':

    cat1 = sys.argv[1]
    cat2 = sys.argv[2]
    title = sys.argv[3]
    clus = sys.argv[4]
    outdir = sys.argv[5]
    class_star = float(sys.argv[6])

    ra = 'ALPHAWIN_J2000'
    dec = 'DELTAWIN_J2000'

    c1 = pyfits.open(cat1)[1].data
    r1 = c1[ra]
    d1 = c1[dec]
    c2 = pyfits.open(cat2)[1].data
    r2 = c2[ra]
    d2 = c2[dec]

    h = htm.HTM()
    maxrad = 1.0/3600.0
    m1,m2,radius = h.match(r1,d1,r2,d2,maxrad,maxmatch=1)
    m1,m2,radius = one_to_one_closest(m1,m2,radius) 

    matched1 = c1[m1]
    matched2 = c2[m2]
    cut = (matched1['CLASS_STAR']>class_star) & (matched2['CLASS_STAR']>class_star)
    matched1 = matched1[cut]
    matched2 = matched2[cut]
    

    r1_final = matched1[ra]
    d1_final = matched1[dec]
    r2_final = matched2[ra]
    d2_final = matched2[dec]


    angdist = np.empty( len(r1_final) )
    for i in range(len(r1_final)):
        orig = ( float(r1_final[i]), float(d1_final[i]) )
        dest = ( float(r2_final[i]), float(d2_final[i]) )
        angdist[i] = slew_angle( orig,dest ) * 3600.0 * 1000.0 # mas


    bmax = np.amax(angdist)
    bins = np.linspace(0,bmax, num=80)
    avg_angs = np.average(angdist)
    std_angs = np.std(angdist)
    plt.figure()
    n,b,p = plt.hist( angdist, bins=bins, normed=1, cumulative=True )
    for i in range(len(b)):
        if n[i] >= 0.9:
            n50 = (b[i+1] + b[i]) / 2.0
            break


    cdir = os.path.join(outdir,'band2band_coadd')
    if not os.path.lexists(cdir):
        os.makedirs(cdir)

    plt.figure()
    bins = np.linspace(0, 2.2*n50, num=40)
    n,b,p = plt.hist( angdist, bins=bins )
    nmax = np.amax(n)
    plt.xlabel( 'Angular offset (mas)' )
    plt.title(title+ ' N=%i'%(len(angdist)))
    plt.text( 2.0/3.0 * b[-1], 3.0/4.0 * nmax, r'$\mu=$' + ' %.1f mas\n'%(avg_angs) + r'$p_{0.9}=$' + ' %.1f mas'%(n50) )
    save = os.path.join( cdir, '%s_hist.pdf' %title)
    plt.savefig( save )


    rdiff = (r1_final - r2_final) * 3600.0 * 1000.0 # mas
    ddiff = (d1_final - d2_final) * 3600.0 * 1000.0 # mas

    angdist = np.empty( len(r1_final) )
    for i in range(len(r1_final)):
        orig = ( float(r1_final[i]), float(d1_final[i]) )
        dest = ( float(r2_final[i]), float(d2_final[i]) )
        angdist[i] = slew_angle( orig,dest ) * 3600.0 * 1000.0 # mas


    rmin = np.amin(r1_final)
    rmax = np.amax(r1_final)
    dmin = np.amin(d1_final)
    dmax = np.amax(d1_final)
    plt.figure()
    plt.xlabel( r'RA ($^{\circ}$)' )
    plt.ylabel( r'DEC ($^{\circ}$)' )
    plt.title(title)
    plt.xlim( [rmin-0.2, rmax+0.2])
    plt.ylim( [dmin-0.2, dmax+0.2])
    qscale = 3000
    quiver( r1_final, d1_final, rdiff, ddiff, pivot='middle', units='width', scale_units='width', scale=qscale, width=0.002 )
    quiver( [rmax],[dmax+0.15], [50], [0], color='b', pivot='middle', units='width', scale_units='width', scale=qscale, width=0.002 )
    plt.text( rmax,dmax+0.12, s='50 mas',color='b',size='small', ha='center', va='center')
    save = os.path.join(cdir, '%s_whisker.pdf' %title)
    plt.savefig( save )

