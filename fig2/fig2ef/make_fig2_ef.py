import sys
sys.path.append("../..")
from bonds_n_coords import coordinates
from bonds_n_coords import lat_bonds
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os
import numpy as np
from matplotlib.ticker import ScalarFormatter

p=3
sc=1000
Nx=[12,9]
Ny = 6
N = [72,54]
U='inf'
SiSj=[]
nh=[0.25,0.33]
center = [33,21] #[int(sys.argv[1]),int(sys.argv[2])]
SiSc=[]
site_num_for_c=[]
fnames=['../../raw_data/file_dmrg_12x6_XC_Uinf_18h_27up_27dn.txt','../../raw_data/file_dmrg_9x6_XC_Uinf_18h_18up_18dn.txt']
bulk_bonds=[]
Rxy=[]
for ff in range(len(fnames)):
    ln = N[ff]*N[ff]
    SiSc.append([])
    site_num_for_c.append([])
    bulk_bonds.append([])
    Rxy.append([])
    cy_bond, op_bond = lat_bonds(Nx[ff],Ny)
    bonds=op_bond
    Rxy[ff] = coordinates(Nx[ff],Ny)
    for i in range(len(bonds)):
        if ((Nx[ff] <= 9 and (bonds[i][0] <= Ny or bonds[i][1] > N[ff]-Ny)) or (Nx[ff] > 9 and (bonds[i][0] <= 2*Ny or bonds[i][1] > N[ff]-2*Ny)) ):
            continue
        bulk_bonds[ff].append(bonds[i])
    command = "grep -A"+str(ln)+" \"i j <Si^z Sj^z>\" "+fnames[ff]+" > tmp"
    os.system(command)
    SS = [[],[],[],[],[]]
    SiSj.append([])
    f=open("tmp",'r')
    for l,line in enumerate(f):
          line=line.strip()
          line=line.split(" ")
          for n in range(len(op_bond)):
               if (l == (op_bond[n][0]-1)*N[ff]+op_bond[n][1]):
                   sitei = int(line[0])
                   sitej = int(line[1])
                   if (Nx[ff] <= 9 and (sitei <= Ny or sitej <= Ny or sitei > N[ff]-Ny or sitej > N[ff]-Ny)):
                         continue
                   if (Nx[ff] > 9 and (sitei <= 2*Ny or sitej <= 2*Ny or sitei > N[ff]-2*Ny or sitej > N[ff]-2*Ny)):
                         continue
                   #print(line[0],line[1])
                   SS[0].append(int(line[0]))
                   SS[1].append(int(line[1]))
                   SS[2].append(float(line[2]))
                   SS[3].append(float(line[3]))
                   SS[4].append(float(line[4]))
                   val = (float(line[2]) + 0.5 * float(line[3]) + 0.5 * float(line[4]))
                   SiSj[ff].append(np.abs(val))
                   #print(line[0],line[1],line[2],line[3],line[4],val)
          
          if l>(center[ff]-1)*N[ff] and l<=center[ff]*N[ff]:
                   sitei = int(line[0])
                   sitej = int(line[1])
                   if (Nx[ff] <= 9 and (sitei <= Ny or sitej <= Ny or sitei > N[ff]-Ny or sitej > N[ff]-Ny)):
                         continue
                   if (Nx[ff] > 9 and (sitei <= 2*Ny or sitej <= 2*Ny or sitei > N[ff]-2*Ny or sitej > N[ff]-2*Ny)):
                         continue
                   vall = (float(line[2]) + 0.5 * float(line[3]) + 0.5 * float(line[4]))
                   if l == (center[ff]-1)*N[ff] + center[ff] :
                       vall=0
                   SiSc[ff].append(vall)
                   site_num_for_c[ff].append(int(line[1]))
                   #print('hhh',line[0],line[1],line[2],line[3],line[4],vall)
          

    f.close()
Nbonds=[]
for i in range(len(fnames)):
      Nbonds.append(len(SiSj[i]))

ticks=[[(0.1**p)*sc,(0.2**p)*sc],[(0.05**p)*sc,(0.1**p)*sc]]
fig,ax = plt.subplots(1,len(fnames),figsize=(15,4.67), gridspec_kw={'width_ratios': [7.6,8.2]})  #,gridspec_kw={'width_ratios': [4, 8]})#, gridspec_kw={'width_ratios': [4, 7, 8]})figsize=(12,6)
plt.subplots_adjust(wspace=0.0)#, hspace=0.6)
norm = [None] * len(fnames)
cmap = [None] * len(fnames)

# Set normalization and colormap
for i in range(len(fnames)):
    norm[i] = mcolors.Normalize(vmin=(0.0000017)*sc, vmax=(0.00224)*sc) # [i]  #min(SiSj)**p max(SiSj)**p
    cmap[i] = cm.Reds #Greys
    for k,bnd in enumerate(bulk_bonds[i]):
                  s1=bnd[0]
                  s2=bnd[1]
                  strength = (SiSj[i][k]**p)*sc
                  color = cmap[i](norm[i](strength))
                  ax[i].plot([Rxy[i][s1][0],Rxy[i][s2][0]],[Rxy[i][s1][1],Rxy[i][s2][1]],c=color,lw=3) #-ss_avg
    ax[i].axis("off")
    ax[i].axis("equal")
    sm = plt.cm.ScalarMappable(cmap=cmap[i], norm=norm[i])
    sm.set_array([])  # Needed for colorbar without imag
cbar = plt.colorbar(sm, ax=ax[i], shrink=0.7, ticks=[0.5,1.0,1.5,2.0], location='right', pad=0)#-0.2)
fs=20
plt.text(8,4.7,'$\\times 10^{-3}$',fontsize=fs)
plt.text(-8.5,4.6,'(e)',fontsize=fs,style='italic')
plt.text(0.8,4.6,'(f)',fontsize=fs,style='italic')
plt.text(-7,-1,"$|\\langle S_{\\chi} \cdot S_{j} \\rangle|$",fontsize=fs)
ax[0].scatter(8,-0.8,s=42,c='r')
ax[1].scatter(1.5,-0.8,s=40,c='b')
plt.text(1.8,-0.95,': Positive',fontsize=fs,style='italic')
plt.text(-1.8,-0.95,': Negative',fontsize=fs,style='italic')
ax[0].text(Rxy[0][center[0]][0]-0.21,Rxy[0][center[0]][1]-0.15,'$\\chi$',fontsize=fs+5,c='black')
ax[1].text(Rxy[1][center[1]][0]-0.22,Rxy[1][center[1]][1]-0.15,'$\\chi$',fontsize=fs+5,c='black')

plt.text(9.4,1.4,"$|\\langle S_{i} \cdot S_{j} \\rangle|^{"+str(int(p))+"}$",fontsize=18.5,rotation=90)
plt.text(-5.2,4.6,'$n_{h}$ = '+str(nh[0]),fontsize=18.5)
plt.text(3.2,4.6,'$n_{h}$ = '+str(nh[1]),fontsize=18.5)
#-----------------------------------------------------------------------------------------------------------------
colorc=[[],[]]
for i in range(2):
    for j in range(len(SiSc[i])):
        if SiSc[i][j]>0:
            colorc[i].append('b')
        else:
            colorc[i].append('r')
for i in range(2):
    for j in range(len(SiSc[i])):
         ax[i].scatter(Rxy[i][site_num_for_c[i][j]][0],Rxy[i][site_num_for_c[i][j]][1],s=np.abs(SiSc[i][j])*5000,c=colorc[i][j],zorder=10)
#-----------------------------------------------------------------------------------------------------------------

#plt.savefig('fig2_ef.pdf',dpi=600,bbox_inches='tight')
plt.show()

