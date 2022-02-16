import os
import sys
import glob
import argparse

import numpy as np
from penquins import Kowalski

# for jd to hjd conversion
from astropy import time, coordinates as coord, units as u

# the database name in Kowalski, this needs to be updated for each update of the database
#DATABASE = "ZTF_sources_20191101"
#DATABASE = "ZTF_sources_20200401"
#DATABASE = "ZTF_sources_20201201"
DATABASE = 'ZTF_sources_20210401'


def JD2HJD(jd,ra,dec,loc='palomar'):
    # convert JD to HJD for give ra and dec, for 
    objectcoords = coord.SkyCoord(ra*u.deg,dec*u.deg, frame='icrs')
    palomar = coord.EarthLocation.of_site(loc)
    times = time.Time(jd, format='jd',scale='utc', location=palomar)

    ltt_helio = times.light_travel_time(objectcoords, 'heliocentric')
    times_heliocentre = times.utc + ltt_helio

    return times_heliocentre.jd



def mag2flux(mag,dmag=[],flux_0 = 3631.0):
    # converts magnitude to flux in Jy asuming AB system
    flux = flux_0 * 10**(-0.4*mag)

    if dmag==[]:
        return flux
    else:
        dflux_p = (flux_0 * 10**(-0.4*(mag-dmag)) - flux)
        dflux_n = (flux_0 * 10**(-0.4*(mag+dmag)) - flux)
        return flux, dflux_p, dflux_n


def get_alertID(ra,dec,kow,radius=1.5):
    # get all alerts for given Ra and Dec

    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": 2,
                "cone_search_unit": "arcsec",
                "radec": {
                    "target": [ra,dec]
                }
            },
            "catalogs": {
                "ZTF_alerts": {
                    "filter": {},
                    "projection": {
                        "_id": 1,
                        "candid": 1,
                        "objectId": 1,
                        "candidate.ra": 1,
                        "candidate.dec": 1
                    }
                }
            }
        },
        "kwargs": {
            "filter_first": False
        }
    }

    #print(q)
    r = kow.query(q)
    #print(r)
    data = r.get('data')

    # 
    key = list(data.keys())[0]
    data = data[key]
    key = list(data.keys())[0]
    data = data[key]

    try:
        return data[0]['objectId'],data[0]['candidate']['ra'],data[0]['candidate']['dec']

    except:
        return None


def get_alertLC(ra,dec,kow,radius=1.5,min_epochs = 1,return_flux=False):
    # get all alerts for given Ra and Dec

    # setup query
    # old
    ##qu = { "query_type": "cone_search", "object_coordinates": { "radec": "[(%.5f,%.5f)]"%(ra,dec), "cone_search_radius": "%.2f"%radius, "cone_search_unit": "arcsec" }, "catalogs": { "ZTF_alerts": { "filter": "{}", "projection": "{'candidate.jd': 1,'candidate.fid': 1, 'candidate.magpsf': 1, 'candidate.sigmapsf': 1, 'candidate.magnr': 1, 'candidate.sigmagnr': 1, 'candidate.distnr': 1, 'candidate.fid': 1, 'candidate.programid': 1, 'candidate.maglim': 1, 'candidate.isdiffpos': 1, 'candidate.ra': 1, 'candidate.dec': 1}" } } }

    ##r = database_query(kow, qu, nquery = 5000)


    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": 2,
                "cone_search_unit": "arcsec",
                "radec": {
                    "target": [ra,dec]
                }
            },
            "catalogs": {
                "ZTF_alerts": {
                    "filter": {},
                    "projection": {
                        "_id": 1,
                        "candid": 1,
                        "objectId": 1,
                        "candidate.jd": 1,
                        "candidate.fid": 1, 
                        "candidate.magpsf": 1,
                        "candidate.sigmapsf": 1,
                        "candidate.magnr": 1, 
                        "candidate.sigmagnr": 1, 
                        "candidate.distnr": 1,
                        "candidate.fid": 1,
                        "candidate.programid": 1,
                        "candidate.maglim": 1, 
                        "candidate.isdiffpos": 1,
                        "candidate.ra": 1,
                        "candidate.dec": 1
                    }
                }
            }
        },
        "kwargs": {
            "filter_first": False
        }
    }

    #print(q)
    r = kow.query(q)
    #print(r)
    data = r.get('data')

    """    
    print(r)
    print(data)
    if not "target" in r:
        print("Query for RA: %.5f, Dec: %.5f failed... returning."%(ra,dec)) 
        return {}
    """

    # 
    key = list(data.keys())[0]
    data = data[key]
    key = list(data.keys())[0]
    data = data[key]

    # storage for outputdata
    magnr,sigmagnr,fid = [],[],[]
    jd, mag, magerr, pos = [], [], [], []
    ra_l, dec_l, fid, pid = [], [], [], []

    # convert the different positive and negative types to 1 and 0
    idp = dict()
    idp['t'] = 1
    idp['1'] = 1
    idp['f'] = 0
    idp['0'] = 0

    # loop over ids to get the data 
    for datlist in data:
        #print(datlist)
        objid = str(datlist["_id"])
        #if not oid is None:
        #    if not objid == str(oid):
        #        continue
        dat = datlist["candidate"]
        jd.append(dat["jd"])
        mag.append(dat["magpsf"])
        magerr.append(dat["sigmapsf"])
        magnr.append(dat["magnr"])
        sigmagnr.append(dat["sigmagnr"])
        pos.append(idp[dat["isdiffpos"]])
        ra_l.append(dat["ra"])
        dec_l.append(dat["dec"])
        fid.append(dat["fid"])
        pid.append(dat["programid"])

    try:
        print(data[0]['objectId'])
    except:
        print('no alerts')



    if np.size(jd)<min_epochs:
        return []

    # convert times to hjd
    jd = np.array(jd)+15./3600/24 # convert from startexp to midexp
    hjd = JD2HJD(jd,ra_l,dec_l)

    # make arrays
    mag = np.array(mag)
    magerr = np.array(magerr)
    magnr = np.array(magnr)
    sigmagnr = np.array(sigmagnr)
    pos = np.array(pos,dtype=float)

    if return_flux:
        # convert mag to flux
        fluxref,fluxref_err,_ = mag2flux(magnr,sigmagnr)
        fluxdiff,fluxdiff_err,_ = mag2flux(mag,magerr)
        flux = fluxref + (-1)**(1.-pos)*fluxdiff
        fluxerr = np.sqrt(fluxdiff_err**2 + fluxref_err**2) # this is conservative, use a minus sign for theoratical errorbar, see Franks document

        # combine everything in a lightcurve    
        lightcurve = np.c_[hjd,flux,fluxerr,fid,pid,ra_l,dec_l,np.ones_like(hjd),np.zeros_like(hjd)]

    else:
        lightcurve = np.c_[hjd,mag,magerr,fid,pid,ra_l,dec_l,np.ones_like(hjd),np.zeros_like(hjd)]
    return lightcurve



def get_PSFlc(ra, dec, kow, radius = 1.5, oid = None, program_ids = [1,2,3], min_epochs = 1,return_flux=False,):

    # setup query
    #qu = { "query_type": "cone_search", "object_coordinates": { "radec": "[(%.5f,%.5f)]"%(ra,dec), "cone_search_radius": "%.2f"%radius, "cone_search_unit": "arcsec" }, "catalogs": { DATABASE : { "filter": "{}", "projection": "{'data.hjd': 1, 'data.mag': 1, 'data.magerr': 1, 'data.fid': 1, 'data.programid': 1, 'data.maglim': 1, 'data.ra': 1, 'data.dec': 1, 'data.catflags': 1}" } } }
    #r = database_query(kow, qu, nquery = 10)

    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": 2,
                "cone_search_unit": "arcsec",
                "radec": {
                    "target": [ra,dec]
                }
            },
            "catalogs": {
                DATABASE: {
                    "filter": {},
                    "projection": {
                        "_id": 1,
                        "data.hjd": 1,
                        "data.fid": 1, 
                        "data.filter": 1, 
                        "data.mag": 1,
                        "data.magerr": 1,
                        "data.ra": 1,
                        "data.dec": 1,
                        "data.programid": 1,
                        "data.catflags": 1}
                }
            }
        },
        "kwargs": {
            "filter_first": False
        }
    }

    r = kow.query(q)
    data = r.get('data')

    """
    print(r)
    print(data)
    if not "target" in r:
        print("Query for RA: %.5f, Dec: %.5f failed... returning."%(ra,dec)) 
        return {}
    """

    # 
    key = list(data.keys())[0]
    data = data[key]
    key = list(data.keys())[0]
    data = data[key]
    #print(data)

    # storage for outputdata
    hjd, mag, magerr = [], [], []
    ra, dec, fid, pid, catflags = [], [], [], [], []

    # loop over ids to get the data 
    for datlist in data:
        objid = str(datlist["_id"])
        _fid = int(str(objid)[7])
        if not oid is None:
            if not objid == str(oid):
                continue
        dat = datlist["data"]

        for dic in dat:
            if not dic["programid"] in program_ids: continue

            hjd.append(dic["hjd"])
            mag.append(dic["mag"])
            magerr.append(dic["magerr"])
            ra.append(dic["ra"])
            dec.append(dic["dec"])
            fid.append(_fid)
            pid.append(dic["programid"])
            catflags.append(dic["catflags"])

    if return_flux:
        flux,fluxerr,_ = mag2flux(np.array(mag),np.array(magerr))
        lightcurve = np.c_[hjd,flux,fluxerr,fid,pid,ra,dec,np.zeros_like(hjd),catflags]

    else:
        lightcurve = np.c_[hjd,mag,magerr,fid,pid,ra,dec,np.zeros_like(hjd),catflags]
 
    lightcurve = lightcurve[np.argsort(lightcurve[:,0])]

    return lightcurve



def database_query(kow, qu, nquery = 5):
    r = {}
    cnt = 0
    while cnt < nquery:
        r = kow.query(query=qu)
        if "result_data" in r:
            break
        cnt = cnt + 1
    return r



def get_fulllightcurve(ra,dec,K,G,radius=1.5,alert_radius=None,refine_pos=True):

    # get the PSF lightcurve in flux units
    lightcurve = get_PSFlc(ra,dec,G,radius,program_ids=[1,2,3],return_flux=True)

    # before getting the alerts, refine the position based on the PSF lightcurve
    if refine_pos and np.size(lightcurve[:,0]>5.):
        med_ra = np.nanmedian(lightcurve[:,5])
        med_dec = np.nanmedian(lightcurve[:,6]) 
        _r = radius if alert_radius is None else alert_radius
        alerts = get_alertLC(med_ra,med_dec,K,_r,return_flux=True)
    else:
        alerts = get_alertLC(ra,dec,K,radius,return_flux=True)

    # combine the data
    alerts = np.atleast_2d(alerts)
    if np.size(alerts)>0:
        lightcurve = np.r_[lightcurve,alerts]

    return lightcurve





def get_Kowalski_catalogs():
    ra,dec = 078.48047465128,+29.91849073615
    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": 5,
                "cone_search_unit": "arcsec",
                "radec": {
                    "target": [ra,dec]
                }
            },
            "catalogs": {
                "IGAPS_DR2": {
                }
            }
        },
    }

    r = kow.query(q)
    data = r.get('data')


    print(r)




