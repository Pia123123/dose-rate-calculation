"""
@author: Pia Kahle
"""

import numpy as np


#Choose the isotope by deleting the "#" in front of the line which should be used and adding a "#" in front of the line which is not needed.

#Isotope="Am"
Isotope="Cm"


#Choose with which geometry the volume should be calculated by deleting the "#" in front of the line which should be used and adding a "#" in front of the lines which are not needed.

#cell="cuboid" #only for activity per volume calculation, dose integration is done using a cylinder
cell="tall cylinder"


if Isotope=="Am":
	
    zk=2.98e-12 #decay constant
    OD=0.9 #OD-value
	
    #stopping power alpha-particles
    aN=1000 #data granularity
    
    #5275keV
    aR1=40.76e-4 #maximum range in cm
    rr1=np.loadtxt("alpha_Energien/alpha5275keV.txt", delimiter="\t", usecols=2) #residual range
    r1=aR1 - rr1 #position = maximum range - residual range
    S1=np.loadtxt("alpha_Energien/alpha5275keV.txt", delimiter="\t", usecols=1) #electronic stopping power in MeV/cm^3
    r1, S1=zip(*sorted(zip(r1, S1))) 
    ar1=np.linspace(0,aR1,aN)
    aS1_interpoliert=np.interp(ar1, r1, S1)
    f1=0.8674 #propability of emission

    #5233keV
    aR2=40.27e-4 #maximum range in cm
    rr2=np.loadtxt("alpha_Energien/alpha5233keV.txt", delimiter="\t", usecols=2) #residual range
    r2=aR2 - rr2 #position = maximum range - residual range
    S2=np.loadtxt("alpha_Energien/alpha5233keV.txt", delimiter="\t", usecols=1) #electronic stopping power in MeV/cm^3
    r2, S2=zip(*sorted(zip(r2, S2)))
    ar2=np.linspace(0,aR2,aN)
    aS2_interpoliert=np.interp(ar2, r2, S2)
    f2=0.1146  #propability of emission

    #5181keV
    aR3=39.66e-4 #maximum range in cm
    rr3=np.loadtxt("alpha_Energien/alpha5181keV.txt", delimiter="\t", usecols=2) #residual range
    r3=aR3 - rr3 #position = maximum range - residual range
    S3=np.loadtxt("alpha_Energien/alpha5181keV.txt", delimiter="\t", usecols=1) #electronic stopping power in MeV/cm^3
    r3, S3=zip(*sorted(zip(r3, S3)))
    ar3=np.linspace(0,aR3,aN)
    aS3_interpoliert=np.interp(ar3, r3, S3)
    f3=0.01383  #propability of emission


if Isotope=="Cm":
	
    zk=6.32e-14
    OD=0.55
	
    #stopping power alpha-particles
    aN=1000 # granularity of the interpolated data
    
    #5078keV
    aR1=40.76e-4 #maximum range in cm
    rr1=np.loadtxt("alpha_Energien/alpha5078keV.txt", delimiter="\t", usecols=2) #residual range
    r1=aR1 - rr1 #position = maximum range - residual range
    S1=np.loadtxt("alpha_Energien/alpha5078keV.txt", delimiter="\t", usecols=1)  #electronic stopping power in MeV/cm^3
    r1, S1=zip(*sorted(zip(r1, S1))) 
    ar1=np.linspace(0,aR1,aN)
    aS1_interpoliert=np.interp(ar1, r1, S1)
    f1=0.75 #propability of emission

    #5035keV
    aR2=40.27e-4 #maximum range in cm
    rr2=np.loadtxt("alpha_Energien/alpha5035keV.txt", delimiter="\t", usecols=2) #residual range
    r2=aR2 - rr2  #position = maximum range - residual range
    S2=np.loadtxt("alpha_Energien/alpha5035keV.txt", delimiter="\t", usecols=1) #electronic stopping power in MeV/cm^3
    r2, S2=zip(*sorted(zip(r2, S2)))
    ar2=np.linspace(0,aR2,aN)
    aS2_interpoliert=np.interp(ar2, r2, S2)
    f2=0.1652 #propability of emission
    

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
NB=OD*1.57e9 #number of cells per ml
VL=10 #volume of the solution in cm^3
Avogadro=6.022e23 #avogadro constant
n_V=50*10**(-9) #concentration of 50nM
N_V=n_V*Avogadro #number of unstable nuclei
AV=zk*N_V*10**-3 #activity per volume in Bq pro cm^3
AV_aufnahme=AV*VL/(NB*VL*VB) #activity per volume in a cell after uptake in Bq/cm^3

#arrays to be filled
ataaussen1=[]
ataaussen2=[]
ataaussen3=[]
atainnen1=[]
atainnen2=[]
atainnen3=[]
atainnen_volleReichweite1=[]
atainnen_volleReichweite2=[]
atainnen_volleReichweite3=[]

N=1000 # granularity of data
nu=np.linspace(0, np.pi/2, N)


for i in range(0, N):
	#distance between center of the cell and cell wall
    if nu[i]<=np.arctan(b/a):
        c=a/np.cos(nu[i]) 
    else:
        c=b/np.cos(np.pi/2-nu[i])
    
    ataaussen1.append(radiusintaussen(ar1, aS1_interpoliert, c))
    ataaussen2.append(radiusintaussen(ar2, aS2_interpoliert, c))
    if Isotope=="Am":
	    ataaussen3.append(radiusintaussen(ar3, aS3_interpoliert, c))
	
    atainnen1.append(radiusintinnen(ar1, aS1_interpoliert, c))  
    atainnen2.append(radiusintinnen(ar2, aS2_interpoliert, c))
    if Isotope=="Am":
	    atainnen3.append(radiusintinnen(ar3, aS3_interpoliert, c))

#contributions to extracellular dose from the different startenergies
aD_aussen1=2*AV*0.5*winkelint(nu, ataaussen1)*f1*(1.602e-19)*3600*1e3*1e6 #in Gy/h
aD_aussen2=2*AV*0.5*winkelint(nu, ataaussen2)*f2*(1.602e-19)*3600*1e3*1e6 #in Gy/h
if Isotope=="Am":
	aD_aussen3=2*AV*0.5*winkelint(nu, ataaussen3)*f3*(1.602e-19)*3600*1e3*1e6 #in Gy/h

#contributions to intracellular dose from the different startenergies
aD_innen1=2*AV_aufnahme*0.5*(winkelint(nu, atainnen1))*f1*(1.602e-19)*3600*1e3*1e6 #in Gy/h
aD_innen2=2*AV_aufnahme*0.5*(winkelint(nu, atainnen2))*f2*(1.602e-19)*3600*1e3*1e6 #in Gy/h

if Isotope=="Am":
	aD_innen3=2*AV_aufnahme*0.5*(winkelint(nu, atainnen3))*f3*(1.602e-19)*3600*1e3*1e6 #in Gy/h
	
	aD_aussen = aD_aussen1+aD_aussen2+aD_aussen3
	aD_innen = aD_innen1+aD_innen2+aD_innen3

if Isotope=="Cm":
	aD_aussen = aD_aussen1+aD_aussen2
	aD_innen = aD_innen1+aD_innen2

print("doserate from solution = ", aD_aussen, "Gy/h")
print("doserate from cells = ", aD_innen, "Gy/h")
