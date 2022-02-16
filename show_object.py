import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import sys
import matplotlib
import argparse
from scipy.optimize import minimize, least_squares
import gatspy
import os

from astropy.time import Time 
from astropy import coordinates as coord 
from astropy import units as u
from scipy.ndimage import gaussian_filter1d

from getPS1 import getcolorim
import astrofunctions as af

# change fontsize
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)


def HJD2BJD(HJD,ra,dec,site='palomar'): 
     """ convert a HJD time (utc) in a BJD time in tdb. Note that this correction  
     is not exact because the time delay is calculated at the HJD time, and not  
     the JD time. The difference is very small however, and should not be a  
     problem for 0.01sec timings.""" 
     target = coord.SkyCoord(ra*u.deg,dec*u.deg, frame='icrs') 
     tsite = coord.EarthLocation.of_site(site) 
     times = Time(HJD, format='jd', 
                       scale='utc', location=tsite) 

     # remove heliocentric correction 
     ltt_helio = times.light_travel_time(target, 'heliocentric') 
     JD = (times.utc - ltt_helio).jd # remove heliocentric correction 

     JDtimes = Time(JD, format='jd', 
                       scale='utc', location=tsite) 

     # do barycentric correction 
     ltt_bary = JDtimes.light_travel_time(target, 'barycentric') 

     BJD = JDtimes.tdb+ltt_bary 

     return BJD.jd 


def weighted_average(y,dy):
    w = dy**-2

    mean = np.sum(y*w)/np.sum(w)
    uncertainty = np.sum(w)**(-0.5)

    return mean,uncertainty


def bin_lc(x,y,dy,p=1,t0=0,bmin=0,bmax=1,nbins=20):
    # make the bins
    bins = np.linspace(bmin,bmax,nbins+1)
    midbins = (bins[:-1]+bins[1:])/2

    # calculate the phases and idx
    ph = (x-t0)/p%1
    idx = np.digitize(ph,bins)

    # calculate the statistics per bin
    meds = np.array([np.median(y[idx==k]) for k in np.arange(1,nbins+1)])
    output = np.array([weighted_average(y[idx==k],dy[idx==k]) 
                for k in np.arange(1,nbins+1)])
    means = output[:,0]
    stds = output[:,1]
    
    return np.c_[midbins,meds,means,stds]



def get_mideclipse(lc,p,t0=0,nbins=1000,sigma=5):

    matplotlib.use('qt5agg') 
    smooth = np.ones(nbins)

    for k in [1,2,3]:
        m = lc[:,3]==k
        if np.sum(m)<10:
            continue    

        #plt.plot((lc[m,0]-t0)/p%1,lc[m,1],'.')
        binned = bin_lc(lc[m,0],lc[m,1],lc[m,2],p,t0,nbins=nbins)

        x = binned[:,0]
        y1 = binned[:,1]
        y2 = binned[:,2]
        dy = binned[:,3]

        # interpolate the y values here
        m = np.isnan(y1)
        y1[m] = np.interp(x[m],x[~m],y1[~m])
        y2[m] = np.interp(x[m],x[~m],y2[~m])
        dy[m] = np.interp(x[m],x[~m],dy[~m])

        # smooth with a gaussian
        s1 = gaussian_filter1d(y1,sigma,mode='wrap')
        #s2 = gaussian_filter1d(y1,k*10,mode='wrap')

        # combine with other band
        smooth *= s1

    # return the minimum of the lc
    print('minphase=',x[np.argmin(smooth)])
    newt0 = t0 + p * x[np.argmin(smooth)]

    print(x[np.argmin(smooth)])

    return newt0

def flux2mag(flux,dflux=[],flux_0 = 3631.0):
    # converts flux (Jy) to magnitudes
    mag = -2.5*np.log10(flux/flux_0)

    if dflux == []:
        return mag

    else:
        dmag_p = -2.5*np.log10((flux-dflux)/flux_0) - mag
        dmag_n = -2.5*np.log10((flux+dflux)/flux_0) - mag

        return mag, dmag_p, dmag_n

def remove_hc(t,midday=0):
    rt = np.floor(t+midday)
    a,b,c,d = np.unique(rt,return_index=True,return_inverse=True,return_counts=True)
    mask = ~np.isin(rt,a[d>15])
    return mask



parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('lcfile',nargs='+')
parser.add_argument('--pmin',default=0.05,type=float)
parser.add_argument('--pmax',default=5.1,type=float)
parser.add_argument('--tmin',default=None,type=float)
parser.add_argument('--tmax',default=None,type=float)
parser.add_argument('--p',default=None,type=float)
parser.add_argument('--t0',default=None,type=float)
parser.add_argument('--detrend',default=None,type=float)
parser.add_argument('--gpr',default=None,type=float)
parser.add_argument('--Nterms',default=1,type=int)
parser.add_argument('--showpower', action='store_true',default=False)
parser.add_argument('--showmag', action='store_true',default=False)
parser.add_argument('--refine', action='store_true',default=False)
parser.add_argument('--noalerts', action='store_true',default=False)
parser.add_argument('--double', action='store_true',default=False)
parser.add_argument('--clean', action='store_true',default=False)
parser.add_argument('--noshow', action='store_true',default=False)
parser.add_argument('--bands',default=[1,2,3],type=float,nargs='+',)
parser.add_argument('--skipdone', action='store_true',default=False)
parser.add_argument('--removehc', action='store_true',default=False)
args = parser.parse_args()

# load HR diagram data
HRdir = "/home/jan/Documents/Projects/Gaia/"
h = np.loadtxt(HRdir+'gaia_h.dat')
xedges = np.loadtxt(HRdir+'gaia_xedges.dat')
yedges = np.loadtxt(HRdir+'gaia_yedges.dat')


# loop
filelist = np.atleast_1d(args.lcfile)
filelist.sort()
for filename in filelist:
    if args.skipdone and os.path.isfile(filename.rstrip('.dat')+'.png'):
        print("skipping %s" %filename)
        continue

    data = np.genfromtxt(filename)
    if np.size(data)<8*10:
        print('No enough data')
        continue

    # data is assumed to be in flux. Col1=time,col2=flux,col3=flux_e,col4=filter,col8=flags
    ra = np.median(data[:,5])
    dec = np.median(data[:,6])

    cim = getcolorim(ra,dec,size=256,filters="grz")    


    data[:,0] = HJD2BJD(data[:,0],ra,dec)
    mag,dmag,_ = flux2mag(data[:,1],data[:,2])

    for k in [1,2,3]:
        m = data[:,3] == k
        med = np.nanmedian(data[m,1])
        data[m,1],data[m,2] = data[m,1]/med,data[m,2]/med

    # remove data
    m = np.isin(data[:,3],args.bands)
    if args.clean:
        m *= (data[:,8] == 0)
        m *= ~((data[:,7] == 1)&(data[:,1]<-0.25)) # remove very negative alerts
    if args.noalerts:
        m *= data[:,7] == 0
    if args.removehc:
        m *= remove_hc(data[:,0])


    if np.sum(m)<20:
        continue


    print("Using %d points" %np.sum(m))

    if args.tmin is not None:
        m *= data[:,0]>args.tmin
    if args.tmax is not None:
        m *= data[:,0]<args.tmax


    data = data[m]
    mag,dmag = mag[m],dmag[m]

    t = data[:,0] - 2400000.5
    y = data[:,1]
    dy = data[:,2]
    fid = data[:,3]

    if args.detrend is not None:
        tfit = np.polyfit(t,y,args.detrend)
        y /= np.poly1d(tfit)(t)

    if args.gpr:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
        for k in args.bands:
            m = data[:,3] == k

            kernel = ConstantKernel()*RBF(10,(5,100)) + WhiteKernel()
            gp = GaussianProcessRegressor(kernel=kernel, alpha=0,
                                  normalize_y=True)
            gp.fit(t[m].reshape(-1, 1),y[m].reshape(-1, 1))
            y_pred, y_std = gp.predict(t[m].reshape(-1,1), return_std=True)

            plt.plot(t[m],y[m])
            plt.plot(t[m],y_pred)
            plt.show()

            y[m] /= y_pred.flatten()


    baseline = np.max(t)-np.min(t)

    if args.p is None:
        from gatspy import periodic
        if args.Nterms == 1:
            model = periodic.LombScargleMultibandFast(fit_period=True)
        elif args.Nterms>1:
            model = periodic.LombScargleMultiband(fit_period=True,Nterms_band=args.Nterms)
        model.optimizer.period_range=(args.pmin, np.min([args.pmax,baseline]))
        model.optimizer.first_pass_coverage=10
        model.fit(t, y, dy, fid)
        p = model.best_period
        print(p)
        t0 = get_mideclipse(np.c_[t,y,dy,fid],p,0)

        if args.showpower:
            fig1 = plt.figure()
            periods = np.linspace(args.pmin, args.pmax, 10**5)
            P = model.periodogram(periods)
            plt.plot(periods*24*60,P)
            plt.axvline(p*24*60,0,1,c='r')
            plt.ylabel('power')
            plt.xlabel('period (min)')
            plt.show()
    else:
        p = float(args.p)
        if args.t0 is None and p>0:
            t0 = get_mideclipse(np.c_[t,y,dy,fid],p,0)
    if args.t0 is not None:
        t0 = float(args.t0)
    if args.double:
        p *= 2

    matplotlib.use('qt5agg')  
    nper =2
    showflagged = False
    #fig, (ax1,ax2) = plt.subplots(2)
    fig = plt.figure(figsize=(14,9))
    spec = fig.add_gridspec(3, 4,
        hspace=0.3,wspace=0.1,left=0.08,right=1-0.02,bottom=0.1,top=1-0.05,)
    ax1 = fig.add_subplot(spec[0, :])
    ax2 = fig.add_subplot(spec[1, :])
    ax3 = fig.add_subplot(spec[2, 0])
    ax4 = fig.add_subplot(spec[2, 2:])

    fig.suptitle('ra=%.4f, dec=%.4f' %(ra,dec))

    # make the model
    if args.refine:
        Np = 4000
        mt = np.r_[np.linspace(t0,t0+p*2,Np),np.linspace(t0,t0+p*2,Np)]
        mf = np.r_[np.ones(Np),2*np.ones(Np)]
        fy = EBmodel_multiband(bestoutput.x,mt,mf,filters=[1,2])

    for n in np.arange(0,2):
        for a,c in zip([1,2,3],['C2','C3','purple']):
            # plot model
            if args.refine:
                ax2.plot((mt[mf==a]-t0)/p,fy[mf==a],ls='-',c=c,lw=0.5,zorder=0)
                #ax3.plot((mt[mf==a]-t0)/p,fy[mf==a],ls='-',c=c,lw=0.5,zorder=0)

            for phot,flagged,marker in zip([0,0,1],[False,True,False],['o','o','x']):
                if flagged and not showflagged:
                    continue
                m = (data[:,3]==a)*(data[:,7]==phot)*((data[:,8]==0)!=flagged)
                if np.sum(m)<1:
                    continue
                mfc="None" if flagged else None
                if p>0:
                    ax2.errorbar((t[m]-t0)/p%1+n,y[m],dy[m],
                        marker=marker,c=c,ls='none',markerfacecolor=mfc,ms=2,lw=0.5)
                if args.showmag:
                    m1 = (dmag < 0.25)
                    ax1.errorbar(t[m*m1],mag[m*m1],dmag[m*m1],
                        marker=marker,c=c,ls='none',markerfacecolor=mfc,ms=2,lw=0.5)
                else:
                    ax1.errorbar(t[m],y[m],dy[m],
                        marker=marker,c=c,ls='none',markerfacecolor=mfc,ms=2,lw=0.5)
                #ax3.errorbar((t[m]-t0)/p%1+n,y[m],dy[m],
                #    marker=marker,c=c,ls='none',markerfacecolor=mfc,ms=2,lw=0.5)

    if args.showmag:
        ax1.set_ylabel('Mag')
        ax1.set_ylim(ax1.get_ylim()[::-1])
    else:
        ax1.set_ylabel('Flux')


    ax2.set_ylabel('Flux')
    ax2.set_xlim(0,nper)

    ax2.set_xlabel(r'Phase (\textbf{p=%f d} / p=%f hrs )' %(p,p*24))

    # download ps1 image and show
    ax3.imshow(cim)
    ax3.set_xticks([])
    ax3.set_xticks([], minor=True)
    ax3.set_yticks([])
    ax3.set_yticks([], minor=True)

    # make HR diagram
    try:
        gaia_dist = af.Xmatch_Gaia_dist([ra,],[dec,],max_distance=2.)
        gaia = af.Xmatch_Gaia_edr3([ra,],[dec,],max_distance=2.)
  
        dist = gaia[0][6] - 5*np.log10(gaia_dist[0][0]) + 5
        dist_b = gaia[0][6] - 5*np.log10(gaia_dist[0][1]) + 5
        dist_B = gaia[0][6] - 5*np.log10(gaia_dist[0][2]) + 5
        
        err = abs(np.array([dist-dist_b,dist-dist_B]))

        
        ax4.errorbar(gaia[0][7]-gaia[0][8],dist,np.array([err]).T,
            fmt='ro',lw=0.5)

        cmin = 10
        cmax = None
        if cmin is not None:
            h[h < cmin] = None
        if cmax is not None:
            h[h > cmax] = None

        f = ax4.pcolormesh(xedges, yedges, h.T)
        ax4.set_xlim(xedges[0], xedges[-1])
        ax4.set_ylim(yedges[0], yedges[-1])

        # fill the rest with scatter (set rasterized=True if saving as vector graphics)
        #ax.scatter(bp_rp, mg, alpha=0.05, s=1, color='k', zorder=0)
        ax4.invert_yaxis()

        # text
        ax4.text(0.99, 0.98, r'plx=$%3.3f \pm %3.3f$' %(gaia[0][0],gaia[0][1]),
             ha='right',
             va='top',
             transform = ax4.transAxes)


    except:
        pass

    ax4.set_ylim(17,-3)
    ax4.set_xlim(-0.9,5.1)

    plt.savefig(filename.rstrip('.dat')+'.png')
    if not args.noshow:
        plt.show()
    plt.close()


