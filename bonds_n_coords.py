import matplotlib.pyplot as plt
import numpy as np


def lat_bonds(Nx,Ny):
    N = Nx*Ny
    bond=[]
    bonds = []
    for i in range(1,N):
        if i%Ny==1 and i<N-Ny:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bonds.append(b)
            bond.append(b)
            b=[i,i+Ny+1]
            bonds.append(b)
            bond.append(b)
            b=[i,i+Ny-1]
            #bonds.append(b)
            bond.append(b)
            b=[i,i+2*Ny-1]
            #bonds.append(b)
            bond.append(b)
   
        elif (i%Ny==0 and i<N):
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)

        elif i>N-Ny:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            if i==N-Ny+1:
                b=[i,i+Ny-1]
                bond.append(b)

        elif i%2==1:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny-1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny+1]
            bond.append(b)
            bonds.append(b)

        elif i%2==0:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)
    return bond, bonds     #bond for cylindrical bc and bonds for open

"""
bond,_ = lat_bonds()

for i  in range(len(bond)):
    print(bond[i])
"""

#finding coordinates
def coordinates(Nx,Ny):
    N = Nx*Ny
    a1 = np.array([1,0])
    a2 = np.array([0.5,np.sqrt(3)/2.0])

    Rxy = {}

    for j in range(Ny,0,-1):
        for i in range(1,Nx+1):
            label = j+(i-1)*Ny
            if j == 6:
                Rxy[label] = (i-1)*a1
            elif j==5:
                Rxy[label] = (i-1)*a1 + a2
            elif j==4:
                Rxy[label] = (i-2)*a1 + 2*a2
            elif j==3:
                Rxy[label] = (i-2)*a1 + 3*a2
            elif j==2:
                Rxy[label] = (i-3)*a1 + 4*a2
            else:
                Rxy[label] = (i-3)*a1 + 5*a2
    return Rxy

#-----------------------------------3sublattice----------------------------------------
def sub(Nx, Ny):
    sublattice = {}
    color = {}
    N = Nx * Ny

    for row in range(1, Ny + 1):
        for col in range(1,Nx+1):
            # Get lattice coordinates (row, col)
            site=row+(col-1)*Ny
            above = site - Ny
            if row % 2 == 1:
                    # Odd-numbe#d62728 row (1-indexed) logic
                    if site == 1:
                          label = 'A'
                          c = '#d62728'
                    elif sublattice.get(above) == 'A':
                          label = 'B'
                          c = '#1f77b4'
                    elif sublattice.get(above) == 'B':
                          label = 'C'
                          c = 'green'
                    else:
                          label = 'A'
                          c = '#d62728'
            else:
                    # Even-numbe#d62728 row logic
                    if site == 2:
                          label = 'B'
                          c = '#1f77b4'
                    elif sublattice.get(above) == 'B':
                          label = 'C'
                          c = 'green'
                    elif sublattice.get(above) == 'C':
                          label = 'A'
                          c = '#d62728'
                    else:
                          label = 'B'
                          c = '#1f77b4'

            sublattice[site] = label
            color[site] = c

    return sublattice, color
#color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
#--------------------------------------------------------------------------------------

'''
Nx = 12
Ny = 6
N = Nx*Ny
bond,bonds = lat_bonds(Nx,Ny)
Rxy = coordinates(Nx,Ny)
sublat, clr= sub(Nx,Ny)
fig = plt.figure(figsize=(20,12))
for k,bnd in enumerate(bonds):
    s1=bnd[0]
    s2=bnd[1]
    
    plt.plot([Rxy[s1][0],Rxy[s2][0]],[Rxy[s1][1],Rxy[s2][1]],lw=2.3,c='k')
    plt.text(Rxy[s1][0]+0.05,Rxy[s1][1],str(s1),fontsize=10)#,weight="bold",color='green')
plt.text(Rxy[N][0]+0.05,Rxy[36][1],str(N),fontsize=10)

#for i in range(1,N+1):
#    plt.scatter(Rxy[i][0],Rxy[i][1],c=clr[i],s=800,zorder=3)

#for i in range(1,N+1):
#    plt.text(Rxy[i][0],Rxy[i][1]+0.1,str(sublat[i]),c=clr[i],fontsize=55)
#-----------------------------------------------temp------------just_for_a_pic-----------------------------
"""
for k,bnd in enumerate(bonds):
    s1=bnd[0]
    s2=bnd[1]
    plt.scatter(Rxy[s1][0],Rxy[s1][1],s=800,c='k')
ll=22

#plt.scatter(Rxy[ll][0],Rxy[ll][1],s=80,c='white')
plt.scatter(Rxy[N][0],Rxy[N][1],s=800,c='k')
"""
plt.text(Rxy[25][0]-1,Rxy[1][1]+0.3,str(Nx)+'$\\times$'+str(Ny),fontsize=55)
plt.axis("off")
plt.axis('equal')
#plt.savefig('XC'+str(Nx)+'x6.pdf',dpi=600,bbox_inches='tight')
print(bond)
print(Rxy)
print(sublat)
print(clr)
plt.show()
'''
