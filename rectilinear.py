#!/usr/bin/env python3
import numpy as np
from scipy.optimize import least_squares 
import matplotlib.pylab as plt
#Fitting funcitons
def tcrit(x,y1,y2):
    return x[0]*(y1-y2)**(1/0.32) + x[1]

def rcrit(x,y1,tc):
    return 2*x[0]*(y1-tc) + 2*x[1]

def fres(x,y,r1,r2,f):
        return f(x,r1,r2)-y
#Functions to calculate VLE with fit data
def vap(t,tc,rc,a1,a2):
    return a2*(t-tc)+rc-0.5*((t-tc)/a1)**(0.32)
def liq(t,tc,rc,a1,a2):
    return a2*(t-tc)+rc+0.5*((t-tc)/a1)**(0.32)

tExt = []
r0 = np.ones(2)
#Read the data in
low = np.loadtxt('vleLow.txt')
hi = np.loadtxt('vleHi.txt')
ind = 5
#fit least squares 
t_lsq = least_squares(fres,r0,args=(low[:ind,0],low[:ind,1],hi[:ind,1],tcrit))
tc = t_lsq.x[1]
r_lsq = least_squares(fres,r0,args=(low[:ind,1]+hi[:ind,1],low[:ind,0],tc,rcrit))
rc = r_lsq.x[1]
print(low)
print(hi)
#print t_lsq,r_lsq

temp = np.linspace(0.03,tc,200)
vr = vap(temp,tc,rc,t_lsq.x[0],r_lsq.x[0])
lr = liq(temp,tc,rc,t_lsq.x[0],r_lsq.x[0])
#plt.subplot(121)
plt.figure(figsize=(5.4,3.5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.errorbar(low[:,1],low[:,0],xerr=low[:,2],fmt='.b')
plt.errorbar(hi[:,1],hi[:,0],xerr=hi[:,2],fmt='.b')
plt.plot(rc,tc,'Dr')
plt.annotate(r"T$_C^*$={:.4f}".format(tc)+"\n"+r"$\rho_C^*$={:.4f}".format(rc),
        xy=(rc,tc), 
		xytext=(22,-50),
		textcoords='offset points', ha='right', va='bottom', 
		bbox=dict(boxstyle='round,pad=0.5', fc="white",alpha=0.5),
		arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0')
        )
plt.plot(vr,temp,'--k',lw=0.8)
plt.plot(lr,temp,'--k',lw=0.8)
plt.xlabel(r'$\rho^*$')
plt.ylabel(r'T$^*$')
plt.xlim(0,0.2)
plt.minorticks_on()
#plt.subplot(122)
#plt.errorbar(low[:,3],low[:,0],xerr=low[:,4],fmt='.')
#plt.errorbar(hi[:,3],hi[:,0],xerr=hi[:,4],fmt='.')
#t_lsq = least_squares(fres,r0,args=(t,cmin[:,2],cmax[:,2],tcrit))
#t_lsqD = least_squares(fres,r0,args=(t,cmin[:,3],cmax[:,3],tcrit))
#tc = t_lsq.x[1]
#tcD = t_lsqD.x[1]
#print(*t_lsq.x)
plt.tight_layout()
#plt.savefig("dgmVLE.png",dpi=600)
plt.show()
