#!/usr/bin/env python
import numpy as np
import sys
from scipy.interpolate import RectBivariateSpline, interp1d

def getdata(arr,key,dict):
    try:
        res=arr[dict[key]]
    except KeyError:
        res=0*arr[0]
    return res

    
if (len(sys.argv)<4):
    print("Format:\n\n  python paintisochrone.py MESA_history_or_isochrone tableATM output\n");
    exit(-1)

# read in the header row of the isochrone file
# check the first 100 rows
with open(sys.argv[1],"r") as f:
    for i in range(100):
        line = f.readline()
        if "mass" in line and len(line)>15 and not ("initial_mass" in line):
            skipnum=i+1
            break
    
history_labels=line.split()

history_label_dict={}
if history_labels[0]=='#':
    for i,l in enumerate(history_labels[1:]):
        la=l.split(':')
        if (len(la)>1):
            history_label_dict[la[-1]]=int(la[0])-1
        else:
            history_label_dict[l]=i
else:
    for i,l in enumerate(history_labels):
        la=l.split(':')
        if (len(la)>1):
            history_label_dict[la[-1]]=int(la[0])-1
        else:
            history_label_dict[l]=i

# print(history_label_dict)
# print(skipnum)

aa = np.loadtxt(sys.argv[1],skiprows=skipnum,unpack=True)
time=getdata(aa,'star_age',history_label_dict)
m1=getdata(aa,'mass',history_label_dict)
logL=getdata(aa,'log_L',history_label_dict)
logLnuc=getdata(aa,'log_Lnuc',history_label_dict)
logTeff=getdata(aa,'log_Teff',history_label_dict)
logR=getdata(aa,'log_R',history_label_dict)
logg=getdata(aa,'log_g',history_label_dict)
logTc=getdata(aa,'log_center_T',history_label_dict)
logsurfz=getdata(aa,'log_surf_cell_z',history_label_dict)
if (logsurfz[0]==0):
    logsurfz=getdata(aa,'log_surf_z',history_label_dict)
if (m1[0]==0):
    m1=getdata(aa,'star_mass',history_label_dict)
    
centerh1=getdata(aa,'center_h1',history_label_dict)
centerhe4=getdata(aa,'center_he4',history_label_dict)
logLHe=getdata(aa,'log_LHe',history_label_dict)
            
# Tvec=10**logTeff

# read in the header row of the spectral file
lines = open(sys.argv[2],"r").readlines()
# skip the first few lines
for l in lines[4:]:
    if l.startswith('!'):
        labels = l.split()
    else:
        break
    
labels = np.array(labels[1:])
nlabels=len(labels)
labels = np.hstack(("star_age","star_mass","log_L",
                    "log_Lnuc","log_Teff","log_R",
                    "log_g","log_surf_z","log_center_T",labels[4:]))

newlabels=[]
olabels=[]
for i,l in enumerate(labels):
    newlabels.append('%d:%s'%(i+1,l))
    olabels.append(l)
labels=newlabels
    
Teff,Logg = np.loadtxt(sys.argv[2],comments='!',unpack=True,usecols=(0,1))
lTeff=np.log10(Teff)
Logguni = np.unique(Logg)
lTeffuni = np.unique(lTeff)
lTeffuni=np.hstack([lTeffuni[0]-0.5, lTeffuni, lTeffuni[-1]+0.5])
xx,yy = np.meshgrid(lTeffuni,Logguni)

Magnitudes = np.loadtxt(sys.argv[2],comments='!',unpack=True,usecols=(range(4,nlabels)))
mag = []
for j in range(len(Magnitudes)):
    mag606 = xx
    for i in range(len(Logguni)):
        #mag606[i] = np.interp(Teffuni,
        #                      Teff[Logg==Logguni[i]],
        #                      Magnitudes[j][Logg==Logguni[i]])
        xf=lTeff[Logg==Logguni[i]]
        yf=Magnitudes[j][Logg==Logguni[i]]
        ii=np.argsort(xf)
        xf=xf[ii]
        yf=yf[ii]
        # mag606[i]=np.interp(lTeffuni,xf,yf)
        ff = interp1d(xf,yf,fill_value="extrapolate")
        mag606[i] = ff(lTeffuni)
    F606 = RectBivariateSpline(Logguni,lTeffuni,mag606)
    F606vec = [float(F606(XX,YY)) for XX,YY in zip(logg,logTeff)]
    F606vec = F606vec - 5*(logR - 7.65) + 5
    mag.append(F606vec)

with open(sys.argv[3],"w") as f:
    f.write(('#'+(' %s'*len(sys.argv))+'\n') % tuple(sys.argv))
    f.write(('#'+(' %s'*len(olabels))+'\n') % tuple(olabels))
    for l in labels:
        f.write('# %s\n'%l)
with open(sys.argv[3],"ab") as f:
    np.savetxt(f,np.transpose(np.vstack((time,m1,logL,
                                               logLnuc,logTeff,logR,
                                               logg,logsurfz,logTc,mag))),
           header=(('%18s    |'+'%19s     |'*6+'%20s     |'+'%20s    |'+'%21s    |'*(len(labels)-9)) %
                   tuple(labels)))
    
