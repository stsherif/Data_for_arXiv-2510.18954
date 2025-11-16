#from bonds_n_coords import lat_bonds
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")   # or "Agg" for non-interactive, "Qt5Agg" if Qt works
import matplotlib.pyplot as plt
import sys
sys.path.append("..")  
from bonds_n_coords import coordinates
import os
import numpy as np
import re
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

folder=Path('.')

print('no inputs')

Nx = 12# int(sys.argv[1])
Ny = 6#int(sys.argv[2])
num = 4#int(sys.argv[3])
fnames=[]
fnames_temp=[]
nh=[]
for f in folder.rglob('*.txt'):
    if f.is_file() and str(Nx)+'x'+str(Ny) in f.name and 'S_of_q' in f.name:
        fnames_temp.append(str(f))

for f in folder.rglob('*.txt'):
    if f.is_file() and str(Nx)+'x'+str(Ny) in f.name and 'occ' in f.name and 'inf' in f.name:
        fnames_temp.append(str(f))

for f in folder.rglob('*.txt'):
    if f.is_file() and str(Nx)+'x'+str(Ny) in f.name and 'occ' in f.name and 'U0' in f.name:
        fnames_temp.append(str(f))
#print(fnames_temp)

for i in range(3):
    Nh_temp=[]
    for j in range(num):
        match = re.search(r'(\d+)h', fnames_temp[j+i*num])
        Nh_temp.append(int(match.group(1)))
    ind = np.argsort(Nh_temp)
    #print(ind)
    for l in ind:
        fnames.append(fnames_temp[l+i*num])

#print(fnames)
#num = int(len(fnames)/3)
N = Nx*Ny

nh = []

#N = Nx*Ny
#print("N",N)
KX=[]
KY=[]
SS=[]
for i in range(num):
    fname = fnames[i]
    KX.append([])
    KY.append([])
    SS.append([])
    f=open(fname,'r')
    for l,line in enumerate(f):
         line=line.strip()
         line=line.split(" ")
         if l == 0:
             Nh = int(line[1])
             nh.append(round(Nh/N,3))
             #max_sq = float(line[2])
         else:
             KX[i].append(float(line[0]))
             KY[i].append(float(line[1]))
             SS[i].append(float(line[4])) #3 for SZ or 4 for Stot
    f.close()
kx=[]
ky=[]
SSF=[]
for i in range(num):
    kx.append([])
    ky.append([])
    SSF.append([])
    kx[i] = np.reshape(KX[i],(2*2*Nx+1,2*2*Ny+1))
    ky[i] = np.reshape(KY[i],(2*2*Nx+1,2*2*Ny+1))
    SSF[i] = np.reshape(SS[i],(2*2*Nx+1,2*2*Ny+1))

KXo=[]
KYo=[]
OCC=[]
for i in range(num):
    fname = fnames[i+num]
    KXo.append([])
    KYo.append([])
    OCC.append([])
    f=open(fname,'r')
    for l,line in enumerate(f):
         line=line.strip()
         line=line.split(" ")
         #if l == 0:
             #Nh = int(line[1])
             #nh.append(round(Nh/N),3)
             #max_sq = float(line[2])
         #else:
         KXo[i].append(float(line[0]))
         KYo[i].append(float(line[1]))
         OCC[i].append(float(line[2]))
    f.close()
kxo=[]
kyo=[]
occ=[]
for i in range(num):
    kxo.append([])
    kyo.append([])
    occ.append([])
    kxo[i] = np.reshape(KXo[i],(2*2*Nx+1,2*2*Ny+1)) #2*Nx+1,2*Ny+1
    kyo[i] = np.reshape(KYo[i],(2*2*Nx+1,2*2*Ny+1))
    occ[i] = np.reshape(OCC[i],(2*2*Nx+1,2*2*Ny+1))    



KXo0=[]
KYo0=[]
OCC0=[]
for i in range(num):
    fname = fnames[i+num*2]
    KXo0.append([])
    KYo0.append([])
    OCC0.append([])
    f=open(fname,'r')
    for l,line in enumerate(f):
         line=line.strip()
         line=line.split(" ")
         KXo0[i].append(float(line[0]))
         KYo0[i].append(float(line[1]))
         OCC0[i].append(float(line[2]))
    f.close()
kxo0=[]
kyo0=[]
occ0=[]
for i in range(num):
    kxo0.append([])
    kyo0.append([])
    occ0.append([])
    kxo0[i] = np.reshape(KXo0[i],(2*2*Nx+1,2*2*Ny+1))
    kyo0[i] = np.reshape(KYo0[i],(2*2*Nx+1,2*2*Ny+1))
    occ0[i] = np.reshape(OCC0[i],(2*2*Nx+1,2*2*Ny+1))

#vx=1.85 #for all
vx=1.53 #for bulk
vn=0.01
vx2=0.97

if Nx==12:
    #vx=1.96 #all
    vx=1.73 #bulk
    #vx=0.73 #sz
    vn=0.03
    vx2=0.99
#---------------------------------------------------------------------------------------------------------
Nx=[6,6,12]
Ny=[6,6,6]

nh2=[]
sq_max=[]
sK=[]
sM=[]
sM2=[]

fnames=['max_sq_sq_at_kpoints_sq_at_mpoints_6x6_bulk.txt','max_sq_sq_at_kpoints_sq_at_mpoints_12x6_bulk.txt']
for i in range(2):
    nh2.append([])
    sq_max.append([])
    sK.append([])
    sM.append([])
    sM2.append([])
    f=open(fnames[i],'r')
    for l,line in enumerate(f):
        if l>0:
            if i==3 and l<2: #i==1
                continue
            line=line.strip()
            line=line.split(" ")
            nh2[i].append(float(line[0]))
            sq_max[i].append(float(line[1]))
            sK[i].append(float(line[2]))
            sM[i].append(float(line[3]))
            sM2[i].append(float(line[4]))
    f.close()
#---------------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PLOT STATIC STRUCTURE FACTOR OVER THE K SPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nrows=3
#print('num',num)
ncols=num
sx=0
sy=0
offx=0.0
if num==3:
    offx=0.65
    sx=5
    sy=5
elif num==4:
    offx=0.68
    sx=6.67
    sy=5
fig= plt.figure(figsize=(19,7)) #17,5.7
gs = gridspec.GridSpec(3, 6, width_ratios=[1.5,1.5,1.5,1.5,1.5,6.5])
ax=[]
plt.rcParams.update({
    #"font.family": "serif",
    #"font.serif": ["Times New Roman"],
    "font.size": 24,
    #'axes.titlesize': 16,     # title size
    'axes.labelsize': 14,     # x/y label size
    'xtick.labelsize': 14,    # x tick size
    'ytick.labelsize': 14,    # y tick size
    'legend.fontsize': 13,    # legend font
    })
for n in range(3):
       for m in range(4):
              ax.append(fig.add_subplot(gs[n, m]))
toshare=fig.add_subplot(gs[0:2, 5])
ax.append(toshare)

ax.append(fig.add_subplot(gs[2, 5],sharex=toshare))              

plt.subplots_adjust(wspace=0.05, hspace=-0.2)
#colors=str(sys.argv[1])
gx=22
gy=12
norm = mcolors.Normalize(vmin=0,vmax=vx)
for i in range(num):
     im1 = ax[i].contourf(kx[i],ky[i],SSF[i],300,cmap=plt.cm.gist_heat, norm=norm)#vmax=vx, vmin=0.0) #vn,extent=[-1, 1, -1, 1] 
     ax1_divider = make_axes_locatable(ax[i])
     im10 = ax[0].contourf(kx[0],ky[0],SSF[0],300,cmap=plt.cm.gist_heat ,norm = norm)# vmax=vx, vmin=0.0)
     ax10_divider = make_axes_locatable(ax[0])
     #for im in im1.collections:
     #       im.set_edgecolor("face")
     ax[i].axis('equal')
     ax[i].set_ylim([-1,1])
     ax[i].set_xlim(-1,1)
     ax[i].set_aspect('equal', adjustable='box')
     #===========================================================================================================
     ax2_divider = make_axes_locatable(ax[i+4])
     im2 = ax[i+4].hexbin(
     kxo0[i].ravel(),
     kyo0[i].ravel(),
     C=occ[i].ravel(),
     gridsize=(gx,gy), #int(N/4),           # number of hexagons across
     cmap=plt.cm.gist_heat,
     vmax=vx2, vmin=0.0)
     im2.set_edgecolor("face")
     ax[i+4].axis('equal')
     ax[i+4].set_aspect('equal', adjustable='box')
     ax[i+4].set_ylim(-1,1)
     ax[i+4].set_xlim(-1,1)
     #==========================================================================================================
     ax3_divider = make_axes_locatable(ax[i+8])
     im3 = ax[i+8].hexbin(
     kxo0[i].ravel(),
     kyo0[i].ravel(),
     C=occ0[i].ravel(),
     gridsize=(gx,gy), #int(N/4),           # number of hexagons across
     cmap=plt.cm.gist_heat,
     vmax=vx2, vmin=0.0)
     im3.set_edgecolor("face")
     ax[i+8].axis('equal')
     ax[i+8].set_aspect('equal', adjustable='box')
     ax[i+8].set_xlim(-1, 1)
     ax[i+8].set_ylim(-1, 1)
#========================color bar=======================================

cax1 = fig.add_axes([0.468, ax[3].get_position().y0+0.04, 0.006, 0.2155])
tick1 = [0.00, round((np.min(SSF[0]) + np.max(SSF[0])) / 2, 2), round(np.max(SSF[0]), 2)-0.01]
cb1 = fig.colorbar(im10, cax=cax1, orientation='vertical', ticks=tick1)
cb1.ax.tick_params(labelsize=14)

cax3 = fig.add_axes([0.468, ax[11].get_position().y0+0.15, 0.006, 0.2155])
tick3 = [0, round((np.min(occ0[num-1]) + np.max(occ0[num-1])) / 2, 2), round(np.max(occ0[num-1]), 2)-0.01]
cb3 = fig.colorbar(im3, cax=cax3, orientation='vertical', ticks=tick3)
cb3.ax.tick_params(labelsize=14)
#=======================================================================
#corners of hexagonal BZ
coord=[[2.0/3,0],[1.0/3,-1.0/np.sqrt(3)],[-1.0/3,-1.0/np.sqrt(3)],[-2.0/3,0],[-1.0/3,1.0/np.sqrt(3)],[1.0/3,1.0/np.sqrt(3)]]
bonds=[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]]

for j,bond in enumerate(bonds):
    site1=bond[0]-1
    site2=bond[1]-1
    for r in range(12):
            ax[r].set_ylim(-1, 1)
            ax[r].tick_params(labelbottom=False, labelleft=False, bottom=False, left=False) 
            ax[r].plot([coord[site1][0],coord[site2][0]],[coord[site1][1],coord[site2][1]],lw=0.9,c='w',linestyle='--')

fig.canvas.draw()
y_line = 0.87
fig_width = fig.get_figwidth()
for i in range(num):
    bbox = ax[i].get_position()
    x_center = (bbox.x0 + bbox.x1) / 2
    fig.text(x_center,y_line+0.015, f'{nh[i]:.3f}',ha='center',va='bottom',fontsize=18)
    fig.text(x_center, y_line, '|', ha='center', va='center', fontsize=18)
fig.lines.append(plt.Line2D([ax[0].get_position().x0+0.02, (ax[3].get_position().x1)-0.005],[y_line, y_line], transform=fig.transFigure, color='black'))
#=============================================================================================
ax_annotate = fig.add_axes([0, 0, 1, 1], zorder=-1)
ax_annotate.axis('off')

# Draw full ray as arrow
start_x = ax[0].get_position().x0+0.011
end_x = ax[3].get_position().x1

ax_annotate.annotate(
    '',
    xy=(end_x+0.0, y_line),        # arrow tip
    xytext=(start_x, y_line),  # arrow tail
    xycoords='figure fraction',
    textcoords='figure fraction',
    arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
)
#=================================================================================================
ax[0].set_ylabel(r'$U/t=\infty$',rotation=90, position=(-1,0.5),fontsize = 20.5)
ax[4].set_ylabel(r'$U/t=\infty$',rotation=90, position=(-1,0.5),fontsize = 20.5)
ax[8].set_ylabel(r'$U/t=0$',rotation=90, position=(-1,0.5),fontsize = 20.5)
plt.text(0.12,y_line+0.015,'$n_{h}$',fontsize = 24)
ax[0].text(-0.12,0.00,r'$\Gamma$',fontsize = 16,color='w')
ax[0].text(-0.12,1/np.sqrt(3),'$M_{1}$',fontsize = 14,color='w')
ax[0].text(1/3,1/np.sqrt(3),r'$K$',fontsize = 14,color='w')
ax[0].text(0.5,1/(2*np.sqrt(3)),r'$M_{2}$',fontsize = 14,color='w')
#=========================================================================
#==========================================================================
#ax[13].set_ylabel('Spin structure factor',fontsize=24)
ax[12].text(-0.073,0,'$S(\\vec{q})$',fontsize=24,rotation=90)
#'Spin structure factor',fontsize=24,rotation=90)
plt.text(0.55,y_line+0,'$n_{h}$',fontsize = 24)

plt.text(0.085,0.8,'(a)',fontsize=22,style='italic')
plt.text(0.518,0.8,'(d)',fontsize=22,style='italic')
plt.text(0.085,0.56,'(b)',fontsize=22,style='italic')
plt.text(0.085,0.325,'(c)',fontsize=22,style='italic')
#---------------------------------------------------------
ms=9
for i in range(2):
    ax[i+12].plot(nh2[i],sq_max[i],label = 'Max $S$($\\vec{q}$)',markersize=ms,marker='o',linestyle='--',markerfacecolor='none')
    ax[i+12].plot(nh2[i],sK[i],label = '$S(K)$',markersize=ms,marker='o',linestyle='--',markerfacecolor='none')
    ax[i+12].plot(nh2[i],sM[i],label = '$S$($M_{1}$)',markersize=ms,marker='o',linestyle='--',markerfacecolor='none')
    ax[i+12].plot(nh2[i],sM2[i],label = '$S$($M_{2}$)',markersize=ms,marker='o',linestyle='--',markerfacecolor='none')
    ax[i+12].xaxis.tick_top()
    ax[i+12].set_xticks([0.1, 0.2, 0.3, 0.4,0.5,0.6])
    #ax[i].yaxis.tick_right()
    ax[i+12].tick_params(axis='both', which='major', labelsize=18)
    ax[12].legend(ncol=4,fontsize=13.5)
    ax[i+12].grid(alpha=0.3)
    if i == 0:
        ax[i+12].tick_params(labeltop=True)
    if i != 0:
        ax[i+12].tick_params(labeltop=False)
#ax[13].set_xticklabels([])
#ax[13].set_yticklabels([])

ax[12].axvspan(0.22, 0.3, color='grey', alpha=0.3)
ax[13].axvspan(0.22, 0.3, color='grey', alpha=0.3, label='highlighted segment')
ax[12].axvspan(0.46, 0.51, color='grey', alpha=0.3)
ax[13].axvspan(0.46, 0.51, color='grey', alpha=0.3, label='highlighted segment')

ax[13].text(0.09,1.37,'   HS\nphase',fontsize=18,color='navy', style='italic')
ax[13].text(0.3,1.37,'Intermediate \n     phase',fontsize=18,color='navy', style='italic')
ax[13].text(0.515,1.37,'Paramagnetic \n     phase',fontsize=18,color='navy', style='italic')


ax[3].text(0.2,0.6,'$S(\\vec{q})$',c='w',fontsize=16)
ax[7].text(0.2,0.6,'$n_{\\sigma}(\\vec{q})$',c='w',fontsize=16)
ax[11].text(0.2,0.6,'$n_{\\sigma}(\\vec{q})$',c='w',fontsize=16)

ax[12].text(0.6,1,str(Nx[0])+'x'+str(Ny[0]), fontsize=20, fontweight='bold', color='black')
#ax[13].text(0.6,2,str(Nx[1])+'x'+str(Ny[1]), fontsize=16, fontweight='bold', color='black')
ax[13].text(0.6,1,str(Nx[2])+'x'+str(Ny[2]), fontsize=20, fontweight='bold', color='black')
#ax[13].text(0.4,2,r'$M_{1} \leftrightarrow M_{2}$',fontsize=12)
pos = ax[12].get_position()
xx= 0.1
ax[12].set_position([pos.x0, pos.y0+0.18 , pos.width, pos.height-0.22])
#pos = ax[13].get_position()
#ax[13].set_position([pos.x0, pos.y0+0.01 , pos.width, pos.height-0.019])
pos = ax[13].get_position()
ax[13].set_position([pos.x0, pos.y0+0.07 , pos.width, pos.height+0.02])

#---------------------------------------------------------
#plt.savefig("fig2_a_d.pdf",dpi=300,bbox_inches='tight',pad_inches=0.0)
plt.show()
