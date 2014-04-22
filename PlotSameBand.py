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

if __name__=='__main__':

    stars_in_coadd_matched = sys.argv[1]
    atleast = int(sys.argv[2])
    clus = sys.argv[3]
    filter = sys.argv[4]
    plotdir = sys.argv[5]
    class_star = float(sys.argv[6])

    if len(sys.argv)>8:
        qscale = float(sys.argv[6])
    else:
        qscale = 3000

    data = pyfits.open(stars_in_coadd_matched)[1].data

    rdist = {}
    ddist = {}
    jrdist = {}
    jddist = {}
    ardist = {}
    addist = {}
    angdist = {}


    cdir = os.path.join( plotdir, 'oneband_SE2coadd' )
    sdir = os.path.join( cdir, '%s'%(filter) )
    idir = os.path.join( sdir,'individual' )
    dirs = [cdir, sdir, idir]
    for dir in dirs:
        if not os.path.exists(dir):
            subprocess.call( ['mkdir', dir] )

    
    cut = (data['SE_INDEX']!=-1)
    ok = (np.sum(cut,axis=1)>=atleast) & (data['CLASS_STAR']>class_star)
    enough = data[ok]
    
    s_ind = enough['SE_INDEX']
    s_ra = enough['SE_RA']
    c_ra = enough['ALPHAWIN_J2000']
    s_dec = enough['SE_DEC']
    c_dec = enough['DELTAWIN_J2000']
   
    angs = []
    for i in range(len(c_ra)):
        orig = ( float(c_ra[i]), float(c_dec[i]) )
        line_cut = ( s_ind[i]!=-1 )

        ras = s_ra[i][line_cut]
        decs = s_dec[i][line_cut]
        for j in range(len(ras)):
            dest = ( float(ras[j]), float(decs[j]) )
            angd = slew_angle( orig,dest ) * 3600.0
            angs.append(angd)

    angs = np.array(angs)
    avg_angs = np.average(angs)
    std_angs = np.std(angs)


    bmax = np.amax(angs)
    bins = np.linspace(0,bmax*0.27*1000.0, num=80)

    plt.figure()
    n,b,p = plt.hist( angs * 1000.0, bins=bins, normed=1, cumulative=True )
    for i in range(len(b)):
        if n[i] >= 0.9:
            n50 = (b[i+1] + b[i]) / 2.0
            break

    bins = np.linspace(0, 2.2*n50, num=40)
    plt.figure()
    n,b,p = plt.hist( angs * 1000.0, bins=bins )
    nmax = np.amax(n)
    plt.xlabel( 'Angular offset (mas)' )
    plt.title('SE offsets from coadd, %s_%s' %(clus,filter) )
    plt.text( 2.0/3.0 * b[-1], 3.0/4.0 * nmax, r'$\mu=$' + ' %.1f mas\n'%(avg_angs*1000.0) + r'$p_{0.9}=$' + ' %.1f mas'%(n50) )
    save = os.path.join( sdir, 'summary.pdf' )
    plt.savefig( save )

    for exp in range(len(data['SE_RA'][0])):
        ra = data['SE_RA']
        dec = data['DELTAWIN_J2000']
        cut = (data['SE_INDEX']!=-1)
        ok = (np.sum(cut,axis=1)>=atleast) & (data['CLASS_STAR']>class_star)
        ma_ra = np.ma.array( ra, mask=-cut )

        avg_ra = data['ALPHAWIN_J2000']

        ra_diff =  ma_ra[:,exp] - avg_ra
        r_diff = ra_diff[-ma_ra.mask[:,exp] & ok ]
        d_diff = dec[-ma_ra.mask[:,exp] & ok ]

        rdist[str(exp)] = r_diff * np.cos(d_diff*np.pi/180.0) * 3600.0 * 1000.0

        jrdist[str(exp)] = ma_ra[:,exp][-ma_ra.mask[:,exp] & ok ].data
        ardist[str(exp)] = avg_ra[-ma_ra.mask[:,exp] & ok ]

        bins = np.linspace(-1.0*0.27*1000.0,1.0*0.27*1000.0, num=60)
        plt.figure()
        plt.hist( rdist[str(exp)], bins=bins )
        plt.xlabel( r'$\cos$(DEC) $\Delta$RA (mas)' )
        plt.title('Exposure %i'%exp )

        save = os.path.join( idir, 'exp%i_RA_hist.pdf'%(exp) )
        plt.savefig( save )

    for exp in range(len(data['SE_DEC'][0])):
        ra = data['SE_DEC']
        cut = (data['SE_INDEX']!=-1)
        ok = (np.sum(cut,axis=1)>=atleast) & (data['CLASS_STAR']>class_star)
        ma_ra = np.ma.array( ra, mask=-cut )

        avg_ra = data['DELTAWIN_J2000']

        ra_diff =  ma_ra[:,exp] - avg_ra
        r_diff = ra_diff[ -ma_ra.mask[:,exp] & ok ]

        ddist[str(exp)] = r_diff * 3600.0 * 1000.0

        jddist[str(exp)] = ma_ra[:,exp][-ma_ra.mask[:,exp] & ok ]
        addist[str(exp)] = avg_ra[-ma_ra.mask[:,exp] & ok ]

        bins = np.linspace(-1.0*0.27*1000.0,1.0*0.27*1000.0, num=60)
        plt.figure()
        plt.hist( ddist[str(exp)], bins=bins )
        plt.xlabel( r'$\Delta$DEC (mas)' )
        plt.title('Exposure %i'%exp )

        save = os.path.join( idir, 'exp%i_DEC_hist.pdf'%(exp) )
        plt.savefig( save )

    for exp in range(len(data['SE_DEC'][0])):
        angdist[str(exp)] = np.empty( len(jrdist[str(exp)]) )
        for j in range(len(jrdist[str(exp)])):
            orig = ( float(ardist[str(exp)][j]), float(addist[str(exp)][j]) )
            dest = ( float(jrdist[str(exp)][j]), float(jddist[str(exp)][j]) )
            angdist[str(exp)][j] = slew_angle( orig,dest ) * 3600.0 * 1000
        bins = np.linspace(-1.0*0.27*1000.0,1.0*0.27*1000.0, num=60)
        plt.figure()
        plt.hist( angdist[str(exp)], bins=bins )
        plt.xlabel( 'Angular offset (mas)' )
        plt.title('Exposure %i'%exp )

        save = os.path.join( idir, 'exp%i_ANG_hist.pdf'%(exp) )
        plt.savefig( save )
    
        plt.figure()
        n,b,p = plt.hist( angdist[str(exp)], bins=bins, normed=1, cumulative=True )
        for i in range(len(b)):
            if n[i] >= 0.9:
                n50 = (b[i+1] + b[i]) / 2.0
                break

        rms = np.std( angdist[str(exp)] / 1000 )
        ag = np.average( angdist[str(exp)] / 1000 )

    for exp in range(len(data['SE_DEC'][0])):
        rmin = np.amin( ardist[str(exp)] )
        rmax = np.amax( ardist[str(exp)] )
        dmin = np.amin( addist[str(exp)] )
        dmax = np.amax( addist[str(exp)] )
    
        plt.figure()
        plt.xlabel( r'RA ($^{\circ}$)' )
        plt.ylabel( r'DEC ($^{\circ}$)' )
        plt.title('Exposure %i'%exp )
       
        plt.xlim( [rmin-0.2, rmax+0.2])
        plt.ylim( [dmin-0.2, dmax+0.2])

        tmpx = ardist[str(exp)][::2]
        tmpy = addist[str(exp)][::2]
        tmpr = rdist[str(exp)][::2]
        tmpd = ddist[str(exp)][::2]
    
        scale = qscale
        quiver( tmpx, tmpy, tmpr, tmpd, pivot='middle', units='width', scale_units='width', scale=scale, width=0.002 )
        quiver( [rmax],[dmax+0.15], [50], [0], color='b', pivot='middle', units='width', width=0.002, scale=scale )
        plt.text( rmax,dmax+0.12, s='50 mas',color='b',size='small', ha='center', va='center')

        save = os.path.join( idir, 'exp%i_whisker.pdf'%(exp) )
        plt.savefig( save )


    for exp in range(len(data['SE_DEC'][0])):
        rmin = np.amin( ardist[str(exp)] )
        rmax = np.amax( ardist[str(exp)] )
        dmin = np.amin( addist[str(exp)] )
        dmax = np.amax( addist[str(exp)] )
        amin = np.amin( angdist[str(exp)] )
        amax = np.amax( angdist[str(exp)] )
    
        plt.figure()
        plt.xlabel( r'RA ($^{\circ}$)' )
        plt.ylabel( r'DEC ($^{\circ}$)' )
        plt.title('Exposure %i'%exp )
       
        plt.xlim( [rmin-0.2, rmax+0.2])
        plt.ylim( [dmin-0.2, dmax+0.2])
        
        plt.scatter( ardist[str(exp)],addist[str(exp)], c=angdist[str(exp)], alpha=0.5, s=80, edgecolor='None', vmin=0, vmax=135)
        cb = plt.colorbar()
        cb.set_label('Angular Offset (mas)')

        save = os.path.join( idir, 'exp%i_off.pdf'%(exp) )
        plt.savefig( save )


    h_file = os.path.join( sdir, 'hist.pdf' )
    w_file = os.path.join( sdir, 'whisker.pdf' )
    a_file = os.path.join( sdir, 'off.pdf' )

    for file in [h_file,w_file,a_file]:
        if os.path.exists(file):
            subprocess.call( ['rm', file] )

    os.system('gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s/*whisker*' %(w_file,idir))
    os.system('gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s/*off*' %(a_file,idir))
    os.system('gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s/*ANG_hist*' %(h_file,idir))

