"""
@author: Pia Kahle
"""

import numpy as np

#choose with which geometry the volume should be calculated

#cell="cuboid" #only for activity per volume calculation, dose integration is done using a cylinder
cell="tall cylinder"

def radiusintaussen (r, S, c): #Integral over the range outside the cell
    atr=[]
    for j in range (0, len(r)-1):
        if r[j]>=c:
            trapez=(r[j+1] - r[j])*(S[j] + S[j+1])/2.0
            atr.append(trapez)
    return sum(atr)
    
def radiusintinnen (r, S, c): #Integral over the range inside the cell
    atr=[]
    for j in range (0, len(r)-1):
        if r[j]<=c:
            trapez=(r[j+1] - r[j])*(S[j] + S[j+1])/2.0
            atr.append(trapez)
    return sum(atr)
    
def winkelint(nu, ta): #Integral over all possible directions in 2D
    atr2=[]
    for k in range(0, N-1):
        trapez2=(nu[k+1] - nu[k])*(np.sin(nu[k])*ta[k] + np.sin(nu[k+1])*ta[k+1])/2.0
        atr2.append(trapez2)
    return sum(atr2)

#parameters
if cell=="cuboid":
	a=0.25e-4 #from center to wall in cm
	b=0.7e-4 #from center to wall in cm
	VB=(2*a*2*a*2*b) #volume of rectangular prism in cm^3
if cell =="tall cylinder":
	a=0.7e-4 #radius in cm
	b=0.25e-4 #half height in cm
	VB=np.pi*b**2*a*2 #volume cylinder in cm^3

OD=0.9 #OD-value
NB=OD*1.57e9 #number of cells per ml
VL=10 #volume of solution in cm^3
T12=2.6234*364*24*60*60 #half-life of Pm-147 in s
zk=np.log(2)/T12 #decay constant in 1/s
Avogadro=6.022e23 #avogadro constant
n_V=50*10**(-9) #concentration of 50nM
N_V=n_V*Avogadro #number of unstable nuclei
AV=zk*N_V*10**-3 #activity per volume in Bq pro cm^3
AV_aufnahme=AV*VL/(NB*VL*VB) #activity per volume in a cell after uptake in Bq/cm^3 

N=1000  #data granularity
nu=np.linspace(0, np.pi/2, N)

#energyspectrum of generated electrons
E=np.loadtxt("spektrum grob gerastert.txt", delimiter="\t", usecols=0)
I=np.loadtxt("spektrum grob gerastert.txt", delimiter="\t", usecols=1)

#arrays to be filled
aeD_innen=[]
aeD_aussen=[]

for i in range (0, len(E)): 
	aN=1000  #granularity of the interpolated data
	rr1=np.loadtxt("elektronen energie reichweite/elektronen energie reichweite " + str(int(E[i]*1000)) + " keV.txt", delimiter="\t", usecols=2) #residual range
	aR1=rr1[len(rr1)-1] #maximum range
	r1=aR1 - rr1 #position = maximum range - residual range
	S1=np.loadtxt("elektronen energie reichweite/elektronen energie reichweite " + str(int(E[i]*1000)) + " keV.txt", delimiter="\t", usecols=1) #electronic stopping power in MeV/cm^3
	r1, S1=zip(*sorted(zip(r1, S1))) 
	ar1=np.linspace(0,aR1,aN)
	aS1_interpoliert=np.interp(ar1, r1, S1)		
	
	#arrays to be filled
	etaaussen1=[]
	etainnen1=[]
	etainnen_volleReichweite1=[]
		
	for j in range(0, N):
		#distance between center of the cell and cell wall
		if nu[j]<=np.arctan(b/a):
			c=a/np.cos(nu[j])
		else:
			c=b/np.cos(np.pi/2-nu[j])        
		
		etaaussen1.append(radiusintaussen(ar1, aS1_interpoliert, c))
		etainnen1.append(radiusintinnen(ar1, aS1_interpoliert, c))        
	
	#contributions to extracellular dose from the different startenergies
	eD_aussen1=2*AV*0.5*winkelint(nu, etaaussen1)*(1.602e-19)*3600*1e3*1e6 #in Gy/h 
	aeD_aussen.append(eD_aussen1*I[i])

	#contributions to intracellular dose from the different startenergies
	eD_innen1=2*AV_aufnahme*0.5*winkelint(nu, etainnen1)*(1.602e-19)*3600*1e3*1e6 #in Gy/h
	aeD_innen.append(eD_innen1*I[i])
	
D_innen=np.sum(aeD_innen)
D_aussen=np.sum(aeD_aussen)

print("doserate from solution = ", D_aussen, "Gy/h")
print("doserate from cells = ", D_innen, "Gy/h")
