import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import math
import matplotlib.pyplot as plt
import sys
import os
import matplotlib.gridspec as gridspec
from pathlib import Path
import re
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.optimize import curve_fit


nh = []
Nx=[6,12]
Ny=[6,6]
N=[36,72]
num_hole_con=[]
num_en=[]
en_dmrg=[]
en_free=[]
en_freen=[]
en_freez=[]
hole_con=[]
ff=['qw_6x6.txt','qw_12x6.txt']
for s in range(len(Nx)):
    num_hole_con.append([])
    num_en.append([])
    g=open('gs_dmrg_energy_'+str(Nx[s])+'x'+str(Ny[s])+'.txt','r')
    for l,line in enumerate(g):
        if l>0:
            line=line.strip()
            line=line.split(" ")
            num_hole_con[s].append(float(line[0]))
            num_en[s].append(float(line[1]))


    fnames=["sorted_energy_U0_"+str(Nx[s])+"x"+str(Ny[s])+"_in_the_first_zone"]

    en_dmrg.append([])
    en_free.append([])
    en_freen.append([])
    en_freez.append([])
    hole_con.append([])
    for q in range(len(fnames)):
        energy_sorted=[]
        name = fnames[q]
        f=open(name,'r')
        for l,line in enumerate(f):
            line=line.strip()
            line=line.split(" ")
            if(l>0):
                 energy_sorted.append(float(line[0]))  
        for p in range(0,N[s],2):
            delta = (N[s]-p)/N[s]
            hole_con[s].append(delta)
            free_energy=0.0
            for i in range(int(p/2)):
                 free_energy += 2*energy_sorted[i]
            en_free[s].append(free_energy)
            en_freen[s].append(free_energy*delta)
            #print("here",delta,N[s]-p,free_energy)
                

    zfiles=[ff[s]]

    for q in range(len(zfiles)):
                   nh.append([])
                   z = []

                   name = zfiles[q]
                   f=open(name,'r')
                   for l,line in enumerate(f):
                         line=line.strip()
                         line=line.split(" ")
                         nh[s].append(float(line[0]))
                         z.append(float(line[1]))
                   for i in range(len(z)):
                         Np = N[s]-round(nh[s][i]*N[s])
                         if Np%2 ==1:
                            Np +=1
                         dn=0.0
                         for j in range(int(Np/2)):
                             dn += 2*energy_sorted[j]
                         #print(nh[s][i]," ",z[i])    
                         en_freez[s].append(dn*z[i])    
    print(z)
#print(nh[1])
#print(en_freez[1])

#-----------------------------------------------------------------------------------------------

folder=Path('.')
print('Provide inputs: no need')
fnames=[]
for f in folder.rglob('*.txt'):
    if f.is_file() and f.name and 'Nq' in f.name:
        fnames.append(str(f))
match = re.search(r'(\d+)x', fnames[0])
Nx_a = int(match.group(1))
match = re.search(r'x(\d+)', fnames[0])
Ny_a = int(match.group(1))
match = re.search(r'(\d+)h', fnames[0])
Nh_a = int(match.group(1))


fnames.append('a_vs_nh_6x6_0qy.txt')
fnames.append('a_vs_nh_10x6_0qy.txt')
fnames.append('a_vs_nh_12x6_0qy.txt')

ky=[]
NN=[]
kx=[]
num_of_xpnts = 1 + 2 * Nx_a
num_of_ypnts = 1 + 2 * Ny_a

f=open(fnames[0],'r')
for l,line in enumerate(f):
      if ((l-(Ny_a+1))%num_of_ypnts == 0):
            line=line.strip()
            line=line.split(" ")
            kx.append(float(line[0]))
            ky.append(float(line[1]))
            NN.append(float(line[2]))
f.close()
#-----------------------------------------------------------------------------------------------
N2=[]
nha=[]
alpha2=[]
Nx=[6,10,12]
Ny=[6,6,6]
for i in range(len(fnames)-1):
    nha.append([])
    alpha2.append([])
    ff=open(fnames[i+1],'r')
    for l,line in enumerate(ff):
            line=line.strip()
            line=line.split(" ")
            if l >= 0:
                nha[i].append(float(line[0]))
                alpha2[i].append(float(line[1]))
    ff.close()
#--------------------------------------------------------------------------------------------------------
nhx=[]
aax=[]
fnames_inset=['a_vs_nh_6x6_0qx.txt','a_vs_nh_10x6_0qx.txt','a_vs_nh_12x6_0qx.txt']
for i in range(len(fnames_inset)):
    nhx.append([])
    aax.append([])
    g=open(fnames_inset[i],'r')
    for l,line in enumerate(g):
            line=line.strip()
            line=line.split(" ")
            if l >= 0:
                nhx[i].append(float(line[0]))
                aax[i].append(float(line[1]))
    g.close()
#---------------------------------------------------------------------------------------------------------
Nx2=[6,12]
Ny2=[6,6]
coord=[-20,-50]
fontz=50
fontz2=100
markz=20
labelz=40
fig= plt.figure(figsize=(30,18))
gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1.5],height_ratios=[2, 2, 8])
ax=[]
#plt.subplots_adjust(hspace=0)
ax.append(fig.add_subplot(gs[0:2, 0]))
tosharex=fig.add_subplot(gs[0, 1])
ax.append(tosharex)
ax.append(fig.add_subplot(gs[1, 1],sharex=tosharex))
ax.append(fig.add_subplot(gs[2, 0]))
ax.append(fig.add_subplot(gs[2, 1],sharex=tosharex))
for i in range(2):
    ax[i+1].text(0.42,coord[i],str(Nx2[i])+'x'+str(Ny2[i]),fontsize=fontz,fontweight='bold')
    ax[i+1].grid(alpha=0.7)

ax[1].scatter(num_hole_con[0],num_en[0],label="$U=\infty$",color='black',s=300)
ax[1].plot(hole_con[0],en_free[0],label="$U=0$",linestyle='dashed',color='red',linewidth=6)
ax[1].plot(hole_con[0],en_freen[0],linestyle='dotted',color='black',linewidth=6)
ax[1].plot(nh[0],en_freez[0],linestyle='dashed',color="black",linewidth=6)

ax[2].scatter(num_hole_con[1],num_en[1],color='black',s=300)
ax[2].plot(hole_con[1],en_free[1],linestyle='dashed',color='red',linewidth=6)
ax[2].plot(hole_con[1],en_freen[1],label=" $n_{h}$",linestyle='dotted',color='black',linewidth=6)
ax[2].plot(nh[1],en_freez[1],label="z",linestyle='dashed',color="black",linewidth=6)
ax[1].tick_params(labelsize=labelz)
ax[2].tick_params(labelsize=labelz)


nh = []
q = []
Nx=[6,10,12]
Ny=[6,6,6]
fnames=['qw_6x6.txt','qw_10x6.txt','qw_12x6.txt']
#fnames=[sys.argv[1],sys.argv[2]]
for i in range(3):
    nh.append([])
    q.append([])
    f=open(fnames[i],'r')
    for l,line in enumerate(f):
        elements = line.strip().split()
        line=line.strip()
        line=line.split(" ")
        nh[i].append(float(line[0]))
        q[i].append(float(line[1]))
    f.close()
#fig,ax=plt.subplots()
alph=[0.4,0.8,0.4]
colorintensity=[1,0.7,1]
for i in range(3):
    ax[0].plot(nh[i],q[i],label = str(Nx[i])+'x'+str(Ny[i]),marker='o',markersize=markz,alpha=colorintensity[i])#+' and U = '+U[i])
ax[0].set_xlabel('$n_{h}$',fontsize=65)
ax[0].set_ylabel('$Z$',fontsize=65)
ax[0].grid(alpha = 0.7)
ax[0].set_yticks(np.arange(0,1,0.2))
ax[0].tick_params(labelsize=labelz)
plt.tight_layout()
ax[0].set_xlim(0,1)
ax[0].set_ylim(0,0.8)
#plt.text(-0.2,0.74,'(b)',fontsize=fontz)


#ax[2].set_xlabel("$n_{h}$",fontsize=28)
plt.text(-0.085,1.02,"Energy",fontsize=55,rotation=90)
#ax[0].set_ylabel("Energy",fontsize=28)
ax[0].legend(loc='lower right',fontsize=45)#,loc='upper left')
ax[1].legend(ncol=1,loc='lower right',fontsize=33)
ax[2].legend(ncol=1,loc='lower right',fontsize=42)
#--------------------------------------------------------------------------
ax[3].plot(kx,NN,label='DMRG',markersize=markz, marker='o',linewidth=6,c='violet')

sub_kx=[]
sub_NN=[]
c=21
sub_kx=kx[c:-c]
sub_NN=NN[c:-c]
print(sub_kx)
kx_fit=[]
for i in range(100):
    #kx_fit.append(-0.15+i*0.003)
    kx_fit.append(-0.3+i*0.003)
kx_fit.append(0.0)
for i in range(101,200):
    #kx_fit.append(-0.15+i*0.003)
    kx_fit.append(-0.3+i*0.003)


#fit function----------------------------- alpha |kx|-------------------------
def func(x, a):
        return a*abs(x)

param1, param_cov1 = curve_fit(func, sub_kx, sub_NN)
#print(param)
fit2 = []
for i in range(len(kx_fit)):
     fit2.append(func(kx_fit[i],param1[0]))
#print(fit2)
ax[3].plot(kx_fit,fit2,label='Linear Fit',linestyle='dashed',linewidth=6,c='navy') #($\\alpha$ = '+str(round(param1[0],3))+')'

#fit function----------------------------- beta kx^2-------------------------
def funcb(x, b):
        return b*x*x

param2, param_cov2 = curve_fit(funcb, sub_kx, sub_NN)
#print(param)
fit3 = []
kx_fit = kx_fit[37:-36]
for i in range(len(kx_fit)):
     fit3.append(funcb(kx_fit[i],param2[0]))
ax[3].plot(kx_fit,fit3,label='Quadratic Fit',linestyle='--',linewidth=6,c='indianred') #($\\beta$ = '+str(round(param2[0],3))+')'

ax[3].set_xlabel('$q_{x}/2\pi$',fontsize = 60)
ax[3].set_ylim(-0.002,0.13)
ax[3].text(-0.42,0.05,'$N(q)$',rotation=90, fontsize = 65)
ax[3].legend(fontsize=40,loc='upper center')
#leg.get_frame().set_alpha(0)
ax[3].grid(alpha=0.7)
ax[3].set_xlim(-0.31,0.31)
ax[3].set_ylim(-0.0001,0.101)
ax[3].text(0.15,0.03,'24x4', fontsize=fontz, fontweight='bold', color='black')
ax[3].text(0.1,0.02,r'$\mathbf{n_{h} = 0.1}$', fontsize=fontz, color='black')
ax[3].set_xticks((-0.3,0.0,0.3))#,fontsize=labelz)
ax[3].set_yticks((0.0,0.05,0.1))#,fontsize=labelz)
#--------------------------------------------------------------------------
for i in range(3):
    ax[4].plot(nha[i],alpha2[i],label = str(Nx[i])+'x'+str(Ny[i]),alpha=colorintensity[i],marker='o',markersize=markz)
ax[4].set_ylim(0,0.75)
ax[4].set_xlim(-0.01,0.94)
ax[4].grid(alpha=0.7)
#ax[4].legend(fontsize=fontz)
ax[4].set_xlabel('$n_{h}$',fontsize = 65)
ax[4].text(-0.08,0.35,'$MW$',rotation=90, fontsize = 65)
#ax[4].set_ylabel('$MW$',fontsize = fontz)
ax[4].set_yticks((0.0,0.2,0.4,0.6))#,fontsize=fontz)
ax[4].set_xticks((0.0,0.2,0.4,0.6,0.8))#,fontsize=fontz)
ax[4].text(0.05,0.15,'      Holes as \n charge carriers', fontsize=60, color='black',zorder=10)
ax[4].text(0.35,0.5,'    Particles as \n charge carriers', fontsize=60, color='black',zorder=10)

#ax[1].set_aspect('equal', adjustable='box')

#-----------------------------------------------------------------------------
pos = ax[1].get_position()
ax[1].set_position([pos.x0+0.01, pos.y0-0.03 , pos.width, pos.height+0.03])
pos = ax[2].get_position()
ax[2].set_position([pos.x0+0.01, pos.y0 , pos.width, pos.height+0.03])
pos = ax[4].get_position()
ax[4].set_position([pos.x0+0.01, pos.y0-0.0 , pos.width, pos.height])
pos = ax[3].get_position()
ax[3].set_position([pos.x0, pos.y0-0.0 , pos.width, pos.height])
ax[1].tick_params(axis='x', labelbottom=False)
ax[2].tick_params(axis='x', labelbottom=False)

ax[4].set_yticks([0.2, 0.4, 0.6])
ax[4].set_yticklabels(['0.2', ' ','0.6'])

ax[3].set_yticks([0.0, 0.05, 0.1])
ax[3].set_yticklabels(['0.0', ' ','0.1'])

plt.text(-0.71,1.28,'(a)',fontsize=55)
plt.text(0.05,1.28,'(b)',fontsize=55)
plt.text(-0.71,0.66,'(c)',fontsize=55)
plt.text(0.05,0.66,'(d)',fontsize=55)
for a in ax:
    a.tick_params(axis='both', labelsize=55)

#____________________________________________________________________________
inax= inset_axes(ax[4], width="30%", height="30%", loc='lower right')
print('here',nhx)
print(aax)
inax.grid(alpha=0.7)
inax.set_ylim(0.0, 0.78)
inax.set_yticks([0.2, 0.4, 0.6])
inax.text(0.22,0.4,'Periodic direction',fontsize=20)
inax.set_yticklabels(['0.2', '0.4','0.6'])
#ccl=['#1f77b4','#2ca02c']
for i in range(3):
    inax.plot(nhx[i],aax[i],marker='o')#,c=ccl[i])#,markersize=markz)
#____________________________________________________________________________
#-----------------------------------------------------------------------------
plt.savefig("fig3.pdf",dpi=300,bbox_inches='tight')
plt.show()


