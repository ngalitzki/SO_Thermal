import sys
import numpy as np
import scipy.constants as PC
import os
import matplotlib.pyplot as plt
import thermal_lib as thermal
import pdb

pi = np.pi
c = 3.e8

f_i=0.
f_f=50000e9
N_s=5000

##300K_stage##
#zotefoam
Te_zote=260.
emit_zote_300=0.9
emit_zote_50=0.4
emit_zote_4=0.01
emit_zote_350=0.01
r_zote=0.490/2

##300K_stage##
#RT-MLI
Te_RT=168.
emit_RT_300=0.9
emit_RT_50=0.4
emit_RT_4=0.01
emit_RT_350=0.01
r_RT=0.500/2
file_RT="rt_mli.txt" #importing file which contains frequency vs transmission
RT=np.loadtxt(file_RT)
(f_RT, T_RT)=(RT[:,0],RT[:,1])
#(f_RT, T_RT)=np.loadtxt(open(file_RT), unpack=True)
f_RT*=1.0e9

##50K stage##
T_edge_50=45.9
#Alumina filter
emit_AF_50=0.9
emit_AF_4=0.05
r_AF_50=0.440/2
Te_AF_50=55
t_AF_50=0.002
file_AF='Alumina_filter_c.txt'
AF=np.loadtxt(file_AF)
(f_AF_50, T_AF_50)=(AF[:,0],AF[:,1])
#(f_AF_50, T_AF_50)=np.loadtxt(open(file_AF), unpack=True)
f_AF_50*=1.0e9

##4K stage##
T_edge_4=3.8

#Field Lens
emit_FL=0.1
r_FL=0.500/2
t_FL = 0.040
kappa_FL=1
Te_FL=5.
file_FL = 'Alumina_lens_c.txt'
FL = np.loadtxt(file_FL)
(f_FL, T_FL)=(FL[:,0],FL[:,1])
f_FL*=1.0e9

#Aperture Lens
emit_AL=0.1
r_AL=0.500/2
t_AL=0.040
kappa_AL=1
Te_AL = 5.0
#file_AL='Alumina_filter.txt'
(f_AL, T_AL)=(f_FL, T_FL)
#(f_AL, T_AL)=np.loadtxt(open(file_AL), unpack_True)
#f_AL*=1.0e9


n_s = 3.106
d_s = 3.16e-3
loss_s = 0.4e-4

n_ar1 = 2.236
d_ar1 = 279e-6
loss_ar1 = 1e-3

n_ar2 = 1.414
d_ar2 = 432e-6
loss_ar2 = 1e-3

#MMF 12icm
emit_12icm=0.01
r_12icm=0.3/2
kappa_12icm=0.1
t_12icm=0.0028

#MMF 8.7icm
emit_87icm=0.01
r_87icm=0.300/2
kappa_87icm=0.1
t_87icm=0.0028

emit_black = 1
r_black = .365/2
Te_black = 6.

#Collimator Lens
emit_CL=0.1
r_CL=0.530/2
t_CL=0.040
kappa_CL=1
Te_CL=5.0
#file_CL='Alumina_filter.txt'
(f_CL, T_CL)=(f_FL, T_FL)

#(f_CL, T_CL)=np.loadtxt(open(file_CL), unpack_True)
#f_CL*=1.0e9

##350mK stage
T_edge_350=0.350

#MMF 6.3icm
emit_63icm=0.1
r_63icm=0.365/2
kappa_63icm=0.1
t_63icm=0.0028
Te_63icm=0.496

##Filter def##
#Alumina_data='Alumina.dat'
#(T_alumina,kappa_alumina) = np.loadtxt( open(Alumina_data), unpack=True)




#MMF 12icm
file_12icm = '14icm.txt'
icm = np.loadtxt(file_12icm)
(f_12icm,T_12icm) = (icm[:,0], icm[:,1])
f_12icm*=1e9
#T_12icm=T_12icm


#for i in range(len(f_12icm)):
 #   print f_12icm[i],T_12icm[i]

#MMF 8.7icm
file_87icm = '12icm.txt'
icm2 = np.loadtxt(file_87icm)
(f_87icm,T_87icm) = (icm2[:,0], icm2[:,1])
f_87icm*=1e9


#MMF 6.3icm
file_63icm = '107icm.txt'
icm3 = np.loadtxt(file_63icm)
(f_63icm,T_63icm) = (icm3[:,0], icm3[:,1])
f_63icm*=1e9

#-----
# define the incident radiation, polarization state, incident angle
incpol=1. #if E is perpendicular to plane-of-incidence (1), if in the plane-of-incidence (-1)
angle_i=0/180.*pi

##calculate##
df=(f_f-f_i)/N_s
f=f_i

#Zotefoam
P_zote=0.
A_zote=pi*r_zote**2

for i in range(N_s):
    f+=df
    P_zote+=A_zote*emit_zote_300*thermal.rho(f,Te_zote)*df
P_tot=P_zote

print '300K stage'
print 'Zotefoam: ',P_zote
print 'Temperature: ',Te_zote
print 'total power: ',P_tot,'\n'

#RTMLI
f=0.
P_zote=0.
P_RT=0.
A_RT=pi*r_RT**2

for i in range(N_s):
    f+=df
    P_zote+=A_zote*emit_zote_300*thermal.rho(f,Te_zote)*thermal.filter(f_RT,T_RT,f)*df
    P_RT+=A_RT*emit_RT_300*thermal.rho(f,Te_RT)*df
    #pdb.set_trace()
P_tot=P_zote+P_RT

print '300K stage'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Temperature: ',Te_RT
print 'total power: ',P_tot,'\n'

#50K stage#
#Alumina filter
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
A_AF_50=pi*r_AF_50**2

for i in range(N_s):
    f+=df
    P_zote+=A_zote*emit_zote_50*thermal.rho(f,Te_zote)*thermal.filter(f_RT,T_RT,f)*thermal.filter(f_AF_50,T_AF_50,f)*df
    P_RT+=A_RT*emit_RT_50*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*df
    P_AF_50+=A_AF_50*emit_AF_50*thermal.rho(f,Te_AF_50)*df
P_tot=P_zote+P_AF_50+P_RT

print '50K stage'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina filter: ',P_AF_50
print 'Temperature: ',Te_AF_50
print 'total power: ',P_tot,'\n'

#4K stage#
#Field Lens
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
P_FL=0.
A_FL=pi*r_FL**2

for i in range(N_s/10):
    f+=df
    P_zote+=A_zote*emit_zote_4*thermal.rho(f,Te_zote)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_RT,T_RT,f)*thermal.filter(f_FL,T_FL,f)*df
    P_RT+=A_RT*emit_RT_4*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*df
    P_AF_50+=A_AF_50*emit_AF_4*thermal.rho(f,Te_AF_50)*thermal.filter(f_FL,T_FL,f)*df
    P_FL+=A_FL*emit_FL*thermal.rho(f,Te_FL)*df
P_tot=P_zote+P_RT+P_AF_50+P_FL

print '4K stage: Field Lens'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina filter: ', P_AF_50
print 'Field Lens: ',P_FL
print 'Temperature: ',Te_FL
print 'total power: ',P_tot,'\n'

#Aperture Lens
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
P_FL=0.
P_AL=0.
A_AL=pi*r_AL**2

for i in range(N_s/10):
    f+=df
    P_zote+=A_zote*emit_zote_4*thermal.rho(f,Te_zote)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_RT,T_RT,f)*df
    P_RT+=A_RT*emit_RT_4*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*df
    P_AF_50+=A_AF_50*emit_AF_4*thermal.rho(f,Te_AF_50)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*df
    P_FL+=A_FL*emit_FL*thermal.rho(f,Te_FL)*thermal.filter(f_AL,T_AL,f)*df
    P_AL+=A_AL*emit_AL*thermal.rho(f,Te_AL)*df
P_tot=P_zote+P_RT+P_AF_50+P_FL+P_AL

print '4K stage: Aperture lens'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina filter: ', P_AF_50
print 'Field Lens: ',P_FL
print 'Aperture Lens: ',P_AL
print 'Temperature: ',Te_AL
print 'total power: ',P_tot,'\n'

#MMF 12icm
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
P_FL=0.
P_AL=0.
P_HWP=0.
P_12icm=0.
A_12icm=pi*r_12icm**2
Te_12icm=4.0

for i in range(N_s/100):
    f+=df
    P_zote+=A_zote*emit_zote_4*thermal.rho(f,Te_zote)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_RT,T_RT,f)*thermal.filter(f_12icm,T_12icm,f)*df
    P_RT+=A_RT*emit_RT_4*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*df
    P_AF_50+=A_AF_50*emit_AF_4*thermal.rho(f,Te_AF_50)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*df
    P_FL+=A_FL*emit_FL*thermal.rho(f,Te_FL)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*df
    P_AL+=A_AL*emit_AL*thermal.rho(f,Te_AL)*thermal.filter(f_12icm,T_12icm,f)*df
    P_12icm+=A_12icm*emit_12icm*thermal.rho(f,Te_12icm)*df
P_tot=P_zote+P_RT+P_AF_50+P_FL+P_AL+P_12icm

print '4K stage: MMF 12icm'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina filter: ',P_AF_50
print 'Field Lens: ',P_FL
print 'Aperture Lens: ',P_AL
print 'MMF 12icm: ',P_12icm
print 'Temperature: ',Te_12icm
print 'total power: ',P_tot,'\n'

#MMF 8.7icm
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
P_FL=0.
P_AL=0.
P_HWP=0.
P_12icm=0.
P_87icm=0.
A_87icm=pi*r_87icm**2
Te_87icm=4.0

for i in range(N_s/100):
    f+=df
    P_zote+=A_zote*emit_zote_4*thermal.rho(f,Te_zote)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_RT,T_RT,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*df
    P_RT+=A_RT*emit_RT_4*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*df
    P_AF_50+=A_AF_50*emit_AF_4*thermal.rho(f,Te_AF_50)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*df
    P_FL+=A_FL*emit_FL*thermal.rho(f,Te_FL)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*df
    P_AL+=A_AL*emit_AL*thermal.rho(f,Te_AL)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*df
    P_12icm+=A_12icm*emit_12icm*thermal.rho(f,Te_12icm)*thermal.filter(f_87icm,T_87icm,f)*df
    P_87icm+=A_87icm*emit_87icm*thermal.rho(f,Te_87icm)*df
P_tot=P_zote+P_RT+P_AF_50+P_FL+P_AL+P_12icm+P_87icm

print '4K stage: MMF 8.7icm'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina filter: ', P_AF_50
print 'Field lens: ',P_FL
print 'Aperture lens :',P_AL
print 'MMF 12icm: ',P_12icm
print 'MMF 8.7icm: ',P_87icm
print 'Temperature: ',Te_87icm
print 'total power: ',P_tot,'\n'


#Collimator Lens
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
P_FL=0.
P_AL=0.
P_HWP=0.
P_12icm=0.
P_87icm=0.
P_CL=0.
A_CL=pi*r_CL**2

for i in range(N_s/100):
    f+=df
    P_zote+=A_zote*emit_zote_4*thermal.rho(f,Te_zote)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_RT,T_RT,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*df
    P_RT+=A_RT*emit_RT_4*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*df
    P_AF_50+=A_AF_50*emit_AF_4*thermal.rho(f,Te_AF_50)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*df
    P_FL+=A_FL*emit_FL*thermal.rho(f,Te_FL)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*df
    P_AL+=A_AL*emit_AL*thermal.rho(f,Te_AL)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*df
    P_12icm+=A_12icm*emit_12icm*thermal.rho(f,Te_12icm)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*df
    P_87icm+=A_87icm*emit_87icm*thermal.rho(f,Te_87icm)*thermal.filter(f_CL,T_CL,f)*df
    P_CL+=A_CL*emit_CL*thermal.rho(f,Te_CL)*df
P_tot=P_zote+P_RT+P_AF_50+P_FL+P_AL+P_12icm+P_87icm+P_CL

print '4K stage: Colimator Lens'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina Filter: ',P_AF_50
print 'Field Lens: ',P_FL
print 'Aperture Lens: ',P_AL
print 'MMF 12icm: ',P_12icm
print 'MMF 8.7icm: ',P_87icm
print 'Colimator Lens: ',P_CL
print 'Temperature: ',Te_CL
print 'total power: ',P_tot,'\n'

#MMF 6.3icm
f=0.
P_zote=0.
P_RT=0.
P_AF_50=0.
P_FL=0.
P_AL=0.
P_HWP=0.
P_12icm=0.
P_87icm=0.
P_CL=0.
P_63icm=0.
A_63icm=pi*r_63icm**2

for i in range(N_s/100):
    f+=df
    P_zote+=A_zote*emit_zote_350*thermal.rho(f,Te_zote)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_RT,T_RT,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    #pdb.set_trace()
    P_RT+=A_RT*emit_RT_350*thermal.rho(f,Te_RT)*thermal.filter(f_AF_50,T_AF_50,f)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    P_AF_50+=A_AF_50*emit_AF_4*thermal.rho(f,Te_AF_50)*thermal.filter(f_FL,T_FL,f)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    P_FL+=A_FL*emit_FL*thermal.rho(f,Te_FL)*thermal.filter(f_AL,T_AL,f)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    P_AL+=A_AL*emit_AL*thermal.rho(f,Te_AL)*thermal.filter(f_12icm,T_12icm,f)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    P_12icm+=A_12icm*emit_12icm*thermal.rho(f,Te_12icm)*thermal.filter(f_87icm,T_87icm,f)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    P_87icm+=A_87icm*emit_87icm*thermal.rho(f,Te_87icm)*thermal.filter(f_CL,T_CL,f)*thermal.filter(f_63icm,T_63icm,f)*df
    P_CL+=A_CL*emit_CL*thermal.rho(f,Te_CL)*thermal.filter(f_63icm,T_63icm,f)*df
    P_63icm+=A_63icm*emit_63icm*thermal.rho(f,Te_87icm)*df*7
P_tot=P_zote+P_RT+P_AF_50+P_FL+P_AL+P_12icm+P_87icm+P_CL+P_63icm

print '350mK stage: MMF6.3icm'
print 'Zotefoam: ',P_zote
print 'RTMLI: ',P_RT
print 'Alumina filter: ', P_AF_50
print 'Field lens: ',P_FL
print 'Aperture lens: ',P_AL
print 'MMF 12icm: ',P_12icm
print 'MMF 8.7icm: ',P_87icm
print 'Colimator lens: ',P_CL
print 'MMF 6.3icm: ',P_63icm
print 'Temperature: ',Te_63icm
print 'total power: ',P_tot,'\n'

f=0.
P_black=0.
A_black=pi*r_black**2

for i in range(N_s/10):
    f+=df
    P_black+=A_black*emit_black*thermal.rho(f,Te_black)*df
P_tot=P_black

print 'Blackbody'
print 'Blackbody: ',P_black
print 'Temperature: ',Te_black
print 'total power: ',P_tot,'\n'

f=0.
P_black=0.
P_12icm = 0.
A_black=pi*r_black**2

for i in range(N_s/50):
    f+=df
    P_black+=A_black*emit_black*thermal.rho(f,Te_black)*thermal.filter(f_63icm,T_63icm,f)*df
    P_63icm+=A_63icm*emit_63icm*thermal.rho(f,Te_63icm)*df
P_tot=P_black+P_63icm

print 'Blackbody filtered'
print 'Blackbody: ',P_black
print 'MMF 6.3icm: ',P_63icm
print 'Temperature: ',Te_63icm
print 'total power: ',P_tot,'\n'

