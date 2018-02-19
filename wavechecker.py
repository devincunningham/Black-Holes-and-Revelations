
# coding: utf-8

# In[1]:

import aplpy
from astropy import units as u, utils
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import gridspec as gridspec, lines as mlines, pyplot as plt
import numpy as np
import pandas as pd
import pyvo as vo


# In[ ]:

scs_optical_results_pd.RA_ICRS


# In[ ]:

for i in range(len(scs_optical_results_pd.columns)):
    print i,":",scs_optical_results_pd.columns[i]


# In[ ]:

scs_radio_results_pd


# In[2]:

def waveseeker(RA=192.491112052,DEC=5.311410068,framesize=4*u.arcmin,frame='icrs',scs_radius = .025):
    
    myLocation = SkyCoord(RA*u.deg, DEC*u.deg, frame = frame)
    scs_radio_query = vo.dal.SCSQuery('https://heasarc.gsfc.nasa.gov/cgi-bin/vo/cone/coneGet.pl?table=nvss&',
                                      pos=(myLocation.ra.deg, myLocation.dec.deg),
                                      radius=scs_radius)
    scs_radio_results = scs_radio_query.execute()
    scs_radio_results_pd = scs_radio_results.votable.to_table().to_pandas()
    scs_radio_results_pd.columns = scs_radio_results.fieldnames
    
    scs_optical_query = vo.dal.SCSQuery('http://wfaudata.roe.ac.uk/sdssdr8-dsa/DirectCone?DSACAT=SDSS_DR8&DSATAB=PhotoObjAll&',
                                        pos=(myLocation.ra.deg, myLocation.dec.deg),
                                        radius=scs_radius)
    scs_optical_results = scs_optical_query.execute()
    scs_optical_results_pd = scs_optical_results.votable.to_table().to_pandas()
    scs_optical_results_pd.columns = scs_optical_results.fieldnames
    
    scs_Xray_query = vo.dal.SCSQuery('http://cda.harvard.edu/cscvo/coneSearch?',
                                     pos=(myLocation.ra.deg, myLocation.dec.deg),
                                     radius=scs_radius)
    scs_Xray_results = scs_Xray_query.execute()
    scs_Xray_results_pd = scs_Xray_results.votable.to_table().to_pandas()
    scs_Xray_results_pd.columns = scs_Xray_results.fieldnames

    return [scs_radio_results_pd, scs_optical_results_pd, scs_Xray_results_pd]


# In[3]:

def waveplotter(RA=192.491112052,DEC=5.311410068,framesize=4*u.arcmin,frame='icrs',scs_radius=.025):
    filename='largegalaxyexample'
    height_ratios = [8]
    width_ratios = [8,8,8]
    wspace, hspace = 0, 0
    
    frame='icrs'
    myLocation = SkyCoord(RA*u.deg, DEC*u.deg, frame = frame)
    sia_urls = ['https://skyview.gsfc.nasa.gov/cgi-bin/vo/sia.pl?survey=nvss&',
                'https://skyview.gsfc.nasa.gov/cgi-bin/vo/sia.pl?survey=sdssdr7&',
                'http://cda.harvard.edu/cscsiap/queryImages?']
    fitsimages = {}
    fitsnames = ['VLA', 'SDSS', 'Chandra']

    # We will use pyVO and astropy to retrieve and open fits images from VLA, SDSS, and Chandra.
    for i, (sia_url, fitsname) in enumerate(zip(sia_urls, fitsnames)):
        # We first use pyVO's Simple Image Access to generate a query.
        sia_query = vo.sia.SIAQuery(sia_url,
                                pos=(myLocation.ra.deg, myLocation.dec.deg),
                                size = framesize, format='image/fits',
                                intersect='covers')
        # We execute this query.
        sia_results = sia_query.execute()
        print "Downloading {} Image...".format(fitsname)
        # We won't always be able to retrieve an image file. When this happens, our executed query returns
        # an unindexable object. 
        try:
            sia_url = sia_results[0].getdataurl()
            fitsimages['{}'.format(fitsname)] = fits.open(utils.data.download_file(sia_url, timeout=300, show_progress=False))
            print "Success!"
        except IndexError:
            fitsimages['{}'.format(fitsname)] = 0
            print "No image currently available for {} :( ".format(fitsname)

        wavedata = waveseeker(RA=RA, DEC=DEC,
                              framesize=framesize,frame=frame,
                              scs_radius=scs_radius)
        scs_radio_results_pd = wavedata[0]
        scs_optical_results_pd = wavedata[1]
        scs_Xray_results_pd = wavedata[2]

    cmap = 'bone'
    names = ['NRAO/VLA Sky Survey - 1.4 GHz', 'Sloan Digital Sky Survey', 'Chandra Source Catalog']
    units = ['Janskies', 'Photons ADU$^{-1}$', 'Photons cm$^{-2}$ s$^{-1}$']

    fig = plt.figure(figsize=(sum(width_ratios) + wspace * (len(width_ratios) - 1),
                              sum(height_ratios) + hspace * (len(height_ratios) - 1)))
    gs = gridspec.GridSpec(len(height_ratios), len(width_ratios),
                           height_ratios=height_ratios, width_ratios=width_ratios)

    for i, image in enumerate([fitsimages['VLA'],
                               fitsimages['SDSS'],
                               fitsimages['Chandra']]):

        ax = fig.add_subplot(gs[i])
        ax.set_title('{}'.format(names[i]),position=[.5, 1.05])

        # Because I'm just throwing aplpy over gridspec subplots, tick marks and tick labels
        # from the subplot axes will sit underneath the aplpy plots. So, let's remove them.
        # While we're at it, let's also remove the splines, which become problematic and
        # annoying when we add colorbars.
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # We only want to plot images when we actually have image data to plot.
        # We have already set empty image files to a zero value. Anything with a nonzero
        # file value will therefore contain an image. Let's plot those.
        if image != 0:
            ax = aplpy.FITSFigure(image,figure=fig,
                                  subplot=list(gs[i].get_position(fig).bounds))
            ax.recenter(myLocation.ra,myLocation.dec,radius=.025)
            ax.show_colorscale(cmap=cmap, stretch='linear', vmid=None)
            if i != 1:
                try:
                    ax.show_contour(levels=10, alpha=1, cmap=cmap)
                except ValueError:
                    print "Unable to produce contours."

            # We have lots of observation data.
            # Let's look at exactly where these observations are in our frame.
            ax.show_markers(myLocation.ra, myLocation.dec, s=1000, marker='+',
                            c='white', linewidth=2, alpha=.5)
            for (ra, dec) in zip(scs_radio_results_pd.ra,scs_radio_results_pd.dec):
                ax.show_markers(ra, dec, s=300, marker='v',
                                facecolor='none', edgecolor='orange', alpha=1,
                                linewidth=2, zorder=20)
            for (ra, dec) in zip(scs_optical_results_pd.ra,scs_optical_results_pd.dec):
                ax.show_markers(ra, dec, s=75, marker='o',
                                facecolor='none', edgecolor='cornsilk', alpha=.75,
                                linewidth=1, zorder=19)
            for (ra, dec) in zip(scs_Xray_results_pd.ra,scs_Xray_results_pd.dec):
                ax.show_markers(ra, dec, s=300, marker='x',
                                facecolor='deeppink', edgecolor='deeppink', alpha=1,
                                linewidth=2, zorder=21)

            # For the life of me I can't change the float precision of these f***ing ticks.
            # If you can figure out a way to do it, please let me know.
            # Anyways, let's set up a colorbar on top of each image.
            # We'll also use the colorbar label to denote units.
            ax.show_colorbar(location='top',
                             pad=0, axis_label_pad=10,
                             axis_label_text='[ {} ]'.format(units[i]))
            ax.colorbar.set_axis_label_font(style='italic')
            ax.colorbar.set_font(style='italic',size='small')

            # A scalebar representing one arcminute.
            ax.add_scalebar(1*u.arcminute, '1 Arcminute', color='white')

            # Let's change the coordinate type of our axes to scalar degrees with a controlled precision.
            ax.set_xaxis_coord_type('scalar')
            ax.set_yaxis_coord_type('scalar')
            ax.tick_labels.set_xformat('%.2f')
            ax.tick_labels.set_yformat('%.2f')
            ax.ticks.hide()

            # We only really need one y axis label since the plots are organized in one row.
            # Let's get rid of the y axis labels for all but the first plot.
            if i != 0:
                ax.axis_labels.hide_y()

        # Smacking down on a legend so everyone knows which markers are which.
        if i == 0:
            VLA_marker = mlines.Line2D([], [], marker='v', linestyle='None',
                                       mfc='none', mec='orange',
                                       markersize=10, label='VLA')
            SDSS_marker = mlines.Line2D([], [], marker='o', linestyle='None',
                                        mfc='none',mec='yellow',
                                        markersize=10, label='SDSS')
            Chandra_marker = mlines.Line2D([], [], marker='x', linestyle='None',
                                           mfc='none', mec='deeppink',
                                           markersize=10, label='Chandra')
            plt.legend(bbox_to_anchor=(0, -4), loc=3,
                       handles=[VLA_marker,SDSS_marker,Chandra_marker],
                       facecolor='whitesmoke', edgecolor='whitesmoke',
                       framealpha=1).set_zorder(102)

        # Every now and then, especially with Chandra (ugh), we will not have image data to plot.
        # We have already set empty image files to a zero value. When this happens,
        # let's let our audience know that no image has been found. It's prettier this way!
        elif image == 0:
            ax = plt.text(0.5, 0.5,'{} Image Not Found'.format(names[i]),
                          horizontalalignment='center',
                          verticalalignment='center',
                          transform = ax.transAxes)
    if filename is not None:
        plt.savefig('{}.png'.format(filename))
    plt.show()


# In[4]:

waveplotter()


# In[ ]:



