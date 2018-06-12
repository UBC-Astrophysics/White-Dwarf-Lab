#!/usr/bin/env python
import numpy as np
import sys
from scipy.interpolate import interp1d
import scipy.integrate as intg

def lum_noconvect(Tc):
    return (Tc/Tc0)**alpha*L0

def dTcdt_noconvect(Tc,t0):
    return(-lum_noconvect(Tc)/cV)

def dTcdt_nofreeze(Tc,t0):
    return(-10**logLfunkTc(-np.log10(Tc))/cV)

if (len(sys.argv)<4):
    print("""
Format:

 python makefakehistory.py history_file start_time fix_time
""")
    exit(-1)

header=''
with open(sys.argv[1],"r") as f:
    for i in range(6):
        header=header+f.readline()

if False:
    logTeffcol=274
    lumcol=275
    logLcol=276
    logRcol=277
    logTccol=281

data=np.genfromtxt(sys.argv[1],skip_header=5,unpack=True,names=True)
simtime=data['star_age']
log_center_T=data['log_center_T']
log_L=data['log_L']

logTcfunk=interp1d(simtime,log_center_T,fill_value='extrapolate')
logLfunk=interp1d(simtime,log_L,fill_value='extrapolate')

logLfunkTc=interp1d(-log_center_T,log_L,fill_value='extrapolate')

starttime=float(sys.argv[2])
fixtime=float(sys.argv[3])
factor=1e-1

alpha=(logLfunk(fixtime*(1+factor))-logLfunk(fixtime*(1-factor)))/(logTcfunk(fixtime*(1+factor))-logTcfunk(fixtime*(1-factor)))
L0=10**logLfunk(fixtime)
Tc0=10**logTcfunk(fixtime)
cV=-10**logLfunk(fixtime)/(10**logTcfunk(fixtime*(1+factor))-10**logTcfunk(fixtime*(1-factor)))*(fixtime*2*factor)

time=np.logspace(7,10,100)
gotime=starttime+time
y0=10**logTcfunk(starttime)

print("# alpha= %g, L0= %g, Tc0= %g cV= %g" % (alpha,L0,Tc0,cV))

Tc_noconvect=intg.odeint(dTcdt_noconvect,y0,gotime)
Tc_nofreeze=intg.odeint(dTcdt_nofreeze,y0,gotime)



print('#   time        tot_time           logL       logLf   logTc  logTcnoc   logLnoc logTcnof  logLnonf dlogTc/dlogt dlogL/dlogT  -L/dTdt')

for t,got,lTnoc,lTnof in zip(time,gotime,np.log10(Tc_noconvect),np.log10(Tc_nofreeze)):
    print('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f' % (t,got,logLfunk(got),
                                                                                                   logLfunkTc(-logTcfunk(got)),
                                                                                                   logTcfunk(got),
                                                                                                   lTnoc,np.log10(lum_noconvect(10**lTnoc)),
                                                                                                   lTnof,logLfunkTc(-lTnof),
                                                                                                   (logTcfunk(got*1.01)-logTcfunk(got*0.99))/(0.02*got)*np.log(10.0)*t,
                           (logLfunk(got*1.01)-logLfunk(got*0.99))/(logTcfunk(got*1.01)-logTcfunk(got*0.99)),
                 -10**logLfunk(got)/(10**logTcfunk(got*1.01)-10**logTcfunk(got*0.99))*(got*0.02)))


gotime=simtime[simtime>starttime]
y0=10**logTcfunk(gotime[0])
logTc_noconvect=np.log10(intg.odeint(dTcdt_noconvect,y0,gotime))
logL_noconvect=np.log10(lum_noconvect(10**logTc_noconvect))
logTc_nofreeze=np.log10(intg.odeint(dTcdt_nofreeze,y0,gotime))
logL_nofreeze=logLfunkTc(-logTc_nofreeze)

for r in zip(gotime,logTc_noconvect,logL_noconvect,logTc_nofreeze,logL_nofreeze):
    print('%g %g %g %g %g' % r)

indexgood=np.arange(len(simtime)).astype(int)
indexgood=indexgood[simtime>starttime]

for i,j in enumerate(indexgood):
    data['log_center_T'][j]=logTc_noconvect[i]
    data['log_L'][j]=logL_noconvect[i]
    data['log_Teff'][j]=(logL_noconvect[i]-2*data['log_R'][j])*0.25+3.76170237
         
with open('%s_noconvect' % sys.argv[1],'wb') as f:
    f.write(header.encode('utf-8'))
    np.savetxt(f,np.transpose(data))

for i,j in enumerate(indexgood):
    data['log_center_T'][j]=logTc_nofreeze[i]
    data['log_L'][j]=logL_nofreeze[i]
    data['log_Teff'][j]=(logL_nofreeze[i]-2*data['log_R'][j])*0.25+3.76170237
         
with open('%s_nofreeze' % sys.argv[1],'wb') as f:
    f.write(header.encode('utf-8'))
    np.savetxt(f,np.transpose(data))
