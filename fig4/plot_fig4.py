import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import matplotlib.gridspec as gridspec


U=[]
sqmax=[]
sk=[]
sisj=[]
a=[]
fnames=[]
Nx=[6,12]
Ny=[6,6]
nh=[0.11,0.11]
N=[24,48]

fnames=['onesitebulk_max_sq_sq_at_kpoints_sq_at_mpoints_6x6_for_site_14.txt','onesitebulk_max_sq_sq_at_kpoints_sq_at_mpoints_12x6_for_site_27.txt','a_vs_U_6x6bulk_0qy.txt','a_vs_U_12x6bulk_0qy.txt','NN_SS_vs_U_per_bonds_bulk_6x6.txt','NN_SS_vs_U_per_bonds_bulk_12x6.txt','NNN_SS_vs_U_per_bonds_bulk_6x6.txt','NNN_SS_vs_U_per_bonds_bulk_12x6.txt','U_d_dbulk_6x6_4h.txt','U_d_dbulk_12x6_8h.txt']

for i in range(2):
    U.append([])
    sqmax.append([])
    sk.append([])
    f=open(fnames[i],'r')
    for l,line in enumerate(f):
        if l>=0:
             elements = line.strip().split()
             line=line.strip()
             line=line.split(" ")
             if line[0] == 'inf':
                  U[i].append(0) #0
                  sqmax[i].append(float(line[1])/N[i])
                  sk[i].append(float(line[2])/N[i])
             else:    
                  U[i].append(1/float(line[0])) #1/
                  sqmax[i].append(float(line[1])/N[i])
                  sk[i].append(float(line[2])/N[i])
                  #sm.append(float(line[3]))
    f.close()

for i in range(2):
    a.append([])
    f=open(fnames[i+2],'r')
    for l,line in enumerate(f):
        if l>0:
             elements = line.strip().split()
             line=line.strip()
             line=line.split(" ")
             if line[0] == 'inf':
                  a[i].append(float(line[1]))
             else:
                  a[i].append(float(line[1]))
    f.close()    
a2=[] 
'''
for i in range(2):
    a2.append([])
    f=open(fnames[i+4],'r')
    for l,line in enumerate(f):
        if l>0:
             elements = line.strip().split()
             line=line.strip()
             line=line.split(" ")
             if line[0] == 'inf':
                  a2[i].append(float(line[1]))
             else:
                  a2[i].append(float(line[1]))
    f.close()
'''
#--------------------------------------------------
for i in range(4):
    sisj.append([])
    f=open(fnames[i+4],'r')
    for l,line in enumerate(f):
        if l>0:
             elements = line.strip().split()
             line=line.strip()
             line=line.split(" ")
             if line[0] == 'inf':
                  #U[i].append(1/float(line[0]))
                  sisj[i].append(float(line[1]))
             else:
                 #U[i].append(1/float(line[0]))
                 sisj[i].append(float(line[1]))
    f.close()
#-----------------------------------------------------

d=[]
d_bulk=[]
for i in range(2):
    d.append([])
    d_bulk.append([])
    f=open(fnames[i+8],'r')
    for l,line in enumerate(f):
        if l>0:
             elements = line.strip().split()
             line=line.strip()
             line=line.split(" ")
             if line[0] == 'inf':
                  d[i].append(float(line[1]))
                  d_bulk[i].append(float(line[2]))
             else:
                 d[i].append(float(line[1]))
                 d_bulk[i].append(float(line[2]))
    f.close()
slop1, intc1=np.polyfit(U[1][-3:],d_bulk[1][-3:],1) #fisrt3
slop2, intc2=np.polyfit(U[1][2:5],d_bulk[1][2:5],1) 

xfit1 = np.linspace(min(U[1]), max(U[1]), 100)
yfit1 = slop1 * xfit1 + intc1

xfit2 = np.linspace(min(U[1]), max(U[1]), 100)
yfit2 = slop2 * xfit2 + intc2


xint = (intc2-intc1)/(slop1-slop2)
Uintercept = round(1/xint,1) #1/

#-----------------------------------------------------
fs=30
ms=11
ls=24
az=1
fig=plt.figure(figsize=(18,13)) #14,13
fig.subplots_adjust(hspace=0.1,wspace=0.4) 
gs = gridspec.GridSpec(2, 2)
ax=[]
toshare1=fig.add_subplot(gs[0, 0])  # Top-left
ax.append(toshare1)
ax.append(fig.add_subplot(gs[0, 1],sharex=toshare1))
toshare2=fig.add_subplot(gs[1, 0])  # Top-left
ax.append(toshare2)
ax.append(fig.add_subplot(gs[1, 1],sharex=toshare2))

ax[2].plot(U[0],sk[0],label=str(Nx[0])+'x'+str(Ny[0]), marker='o', linestyle='-',markersize=ms,alpha=az)#,color='#1f77b4'
ax[2].plot(U[1],sk[1], label=str(Nx[1])+'x'+str(Ny[1]),marker='o', linestyle='-',markersize=ms,alpha=az)
#ax[2].axvline(x=0.01736, color='black', linestyle='--', linewidth=1)
ax[2].text(-0.015,0.025,'      Hole \n Dominated',fontsize=fs-2,color='black', style='italic')#,fontweight='bold')#,rotation=90) #-0.006,1.45 #0.038 #-0.018
ax[2].text(0.03,0.025,'     Spin \n Dominated',fontsize=fs-2,color='black', style='italic')#,fontweight='bold')#0.03,1.45,#0.038
#ax[2].set_ylabel('$S(K)/N$',fontsize=fs)
ax[2].text(-0.038,0.035,'$S(K)/N$',fontsize=fs+2,color='black',rotation=90)
#ax[0].legend(fontsize=fs-2,ncol=2,loc='upper right')

ax[2].tick_params(axis='x', labelsize=ls)
ax[2].tick_params(axis='y', labelsize=ls)
ax[0].tick_params(axis='x', labelsize=ls)
ax[0].tick_params(axis='y', labelsize=ls)
#ax[2].grid(alpha=0.3)
#ax[0].grid(alpha=0.3)
cclor=['#1f77b4', '#ff7f0e']
#-----------------------------------------------------------------------------------
for i in range(2):
    ax[0].plot(U[i],sisj[i], marker='o', linestyle='-',label=str(Nx[i])+'x'+str(Ny[i]),markersize=ms,alpha=az)#,color='#2ca02c')
    ax[1].plot(U[i],sisj[i+2], marker='o', linestyle='-',label=str(Nx[i])+'x'+str(Ny[i]),markersize=ms,alpha=az)
    #ax[i+2].set_xlabel('$t/U$',fontsize=fs)
    #ax[i].grid(alpha=0.3)

#ax[0].text(-0.02,-0.125,'$NN SS$',fontsize=fs,rotation=90)
ax[0].set_ylabel('$\\frac{1}{N} \\sum_{\\langle i,j\\rangle} \\langle S_{i} \cdot S_{i}\\rangle $',fontsize=fs)
#ax[1].text(-0.02,0.066,'$\\frac{1}{N} \\sum_{\\langle \\langle i,j \\rangle \\rangle} \\langle S_{i} \cdot S_{i}\\rangle $',fontsize=fs,rotation=90)
ax[1].set_ylabel('$\\frac{1}{N} \\sum_{\\langle\\langle i,j\\rangle\\rangle} \\langle S_{i} \cdot S_{i}\\rangle $',fontsize=fs)

ax[0].legend(fontsize=fs-2,ncol=2,loc='upper left')
ax[0].tick_params(axis='x', labelsize=ls)
ax[0].tick_params(axis='y', labelsize=ls)
ax[1].tick_params(axis='x', labelsize=ls)
ax[1].tick_params(axis='y', labelsize=ls)
#ax[2].text(-0.028,-0.103,'(c)',fontsize=fs)
ax[0].set_ylim(-0.132,-0.103)

#ax[3].text(-0.023,0.104,'(d)',fontsize=fs)
ax[1].set_ylim(0.04,0.101)

#ax[0].text(0.02,-0.115,'$\\frac{1}{N} \\sum_{\\langle i,j \\rangle \in NN} \\langle S_{i} \cdot S_{i}\\rangle $',fontsize=fs)
#ax[1].text(0.02,0.06,'$\\frac{1}{N} \\sum_{\\langle i,j \\rangle \in NNN} \\langle S_{i} \cdot S_{i}\\rangle $',fontsize=fs)
#-----------------------------------------------------------------------------------
for i in range(2):
    #ax.plot(U[i],d[i], marker='o', linestyle='-',label=str(Nx[i])+'x'+str(Ny[i])+' nh = '+str(nh[i]))#,color='#2ca02c')
    ax[3].plot(U[i],d_bulk[i], marker='o', linestyle='-',label='Bulk '+str(Nx[i])+'x'+str(Ny[i])+' nh = '+str(nh[i]),markersize=ms,alpha=az)#,color='#2ca02c')

ax[3].plot(xfit1,yfit1,alpha=az,color='black', linestyle='--')#, marker='o', linestyle='-',label='Bulk '+str(Nx[i])+'x'+str(Ny[i])+' nh = '+str(nh[i]))
ax[3].plot(xfit2,yfit2,alpha=az,color='black', linestyle='--')#, marker='o', linestyle='-',label='Bulk '+str(Nx[i])+'x'+str(Ny[i])+' nh = '+str(nh[i]))
x=[xint,xint]
y=[0,0.04]
ax[3].tick_params(axis='x', labelsize=ls)
ax[3].tick_params(axis='y', labelsize=ls)
ax[3].text(xint-0.018,0.015,"$U_{c}/t \\simeq 50 $",fontsize=fs+2, style='italic')#str(Uintercept)

ax[3].set_ylim(-0.005,0.041)
#ax[4].text(-0.015,0.043,'(e)',fontsize=fs)

ax[2].set_xlabel('$t/U$',fontsize=fs)
ax[3].set_xlabel('$t/U$',fontsize=fs)
#ax[3].text(-0.02,0.02,'$D$',fontsize=fs,rotation=90)
ax[3].set_ylabel('$\\frac{1}{N} \\sum_{i} \\langle n_{i\\uparrow}n_{i\\downarrow}\\rangle$',fontsize=fs)
#ax[3].text(0.03,0.035,'$D=\\frac{1}{N} \\sum_{i} \\langle n_{i\\uparrow}n_{i\\downarrow}\\rangle $',fontsize=fs)
#ax[3].grid(alpha=0.3)
#-----------------------------------------------------------------------------------
'''
pos = ax[1].get_position()
ax[1].set_position([pos.x0, pos.y0 , pos.width, pos.height])
pos = ax[2].get_position()
ax[2].set_position([pos.x0, pos.y0+0.03 , pos.width, pos.height])
pos = ax[3].get_position()
ax[3].set_position([pos.x0, pos.y0+0.03 , pos.width, pos.height])
'''
ax[0].tick_params(axis='x', labelbottom=False)
ax[1].tick_params(axis='x', labelbottom=False)

#ax[2].set_yticks([0.02,0.04,0.06,0.08,0.09])
#ax[2].set_yticklabels(['0.02','0.04','0.06','0.08',''])
ax[0].set_yticks([-0.13, -0.12, -0.11])
ax[0].set_yticklabels(['-0.13','-0.12', '-0.11'])
#ax[3].set_yticks([0.00, 0.02, 0.04])
#ax[3].set_yticklabels(['0.00',' ', '0.04'])
#ax[1].set_yticks([0.04, 0.06, 0.08,0.1])
#ax[1].set_yticklabels(['0.04','0.06','', '0.1'])
ax[2].set_xticks([0, 0.02,0.04,0.06,0.08,0.1])
#ax[2].set_yticklabels(['-0.13', ' ', '-0.11'])
'''
ax[0].axvline(x=0.01736, color='black', linestyle=':', linewidth=1)
ax[2].axvline(x=0.01736, color='black', linestyle=':', linewidth=1)
ax[3].axvline(x=0.01736, color='black', linestyle=':', linewidth=1)
ax[1].axvline(x=0.01736, color='black', linestyle=':', linewidth=1)
'''
ax[0].axvspan(0.015, 0.025, color='grey', alpha=0.3)
ax[1].axvspan(0.015, 0.025, color='grey', alpha=0.3)
ax[2].axvspan(0.015, 0.025, color='grey', alpha=0.3)
ax[3].axvspan(0.015, 0.025, color='grey', alpha=0.3)

#plt.axis('off')
plt.text(-0.192,0.087,'(a)',fontsize=fs, style='italic')
plt.text(-0.04,0.087,'(b)',fontsize=fs, style='italic')
plt.text(-0.192,0.036,'(c)',fontsize=fs, style='italic')
plt.text(-0.04,0.036,'(d)',fontsize=fs, style='italic')
plt.savefig('fig4.pdf',dpi=300,bbox_inches='tight')
#plt.show()

