#!/home/jan/anaconda3/bin/python
import json
import argparse
import os
import numpy as np
from penquins import Kowalski
import ZTFlcquery



# pass commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename",help='name of file with coordinate')
parser.add_argument("--delimiter",type=str,default=None,help='delimiter of coord file')
parser.add_argument('--cols', nargs=2, default=[0,1])
parser.add_argument('--skiprows', type=int,default=0)
parser.add_argument("--radius","-r",type=float,default=1.5,help='match radius in arcsec')    
parser.add_argument("--refine_pos",action="store_true", default=False,help='refine the position from LC before getting alerts') 
parser.add_argument("--alert_radius",type=float,default=None,help='match radius for alerts; if None, --radius will be used')    
parser.add_argument("--verbose",'-v',action="store_true", default=False,help='verbose') 
parser.add_argument("--savedir",help='directory of where to save files',type=str,default='./')
parser.add_argument('-o',action="store_true", default=False,help='skip existing files')

args = parser.parse_args()

# load the coordinates
coords=np.loadtxt(args.filename,usecols=np.array(args.cols,dtype=int),skiprows=args.skiprows,delimiter=args.delimiter)
coords=np.atleast_2d(coords) # failsafe for the case of only 1 line in file

# setup kowalski interface
#K = Kowalski(username='jroestel',password='Erg$terkWachtwoord')
#G = Kowalski(username='jroestel',password='Erg$terkWachtwoord',host='gloria.caltech.edu')

with open('/home/jan/mysecrets/secrets.json', 'r') as f:
    secrets = json.load(f)

K = Kowalski(**secrets['kowalski'], verbose=False)
G = Kowalski(**secrets['gloria'], verbose=False)


# load the coordinatesa
#coords=np.loadtxt(args.filename,usecols=np.array(args.cols),skiprows=0,delimiter=args.delimiter)
#coords=np.atleast_2d(coords) # failsafe for the case of only 1 line in file

# got and get LCs
for n,(ra,dec) in enumerate(coords):
    if args.verbose:
        print(n,ra,dec)

    outputfilename = '%s/lc_%08.4f_%07.4f.dat' %(args.savedir,ra,dec)

    if os.path.exists(outputfilename) and args.o:
            if args.verbose:
                print('Skipping coords %08.4f %07.4f' %(ra,dec))
            continue

    # get the lightcurve
    lc = ZTFlcquery.get_fulllightcurve(ra,dec,K,G,
        args.radius,args.alert_radius,args.refine_pos)

    header = "HJD \t flux \t error \t fid \t pid \t ra \t dec \t alertphot \t flags "
    np.savetxt(outputfilename,lc,header=header,fmt='%12.12g')


