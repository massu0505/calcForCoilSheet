

#%%
from cmath import tan
from http.client import NETWORK_AUTHENTICATION_REQUIRED
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from plotly.graph_objs import Scatter3d
# from plotly.offline import init_notebook_mode,iplot
# init_notebook_mode()
# from sklearn import linear_model
# from scipy.optimize import curve_fit
import os
import csv
import skrf as rf

#%%
f=6.78*10**6
w=2*np.pi*f
Q1=140
Q2=231
# ke=np.linspace(0.01, 0.4, num=40)
ke=0.2
# l_FromSourceToCp=np.linspace(0, 1.0, num=21)
l_FromSourceToCp=0
Sh_lambda=0.53
Sh_Rsx=120
# Sh_Rsx=np.linspace(10, 150, num=15)
Sh_C=1/(Sh_Rsx*Sh_lambda*f)
Sh_L=Sh_Rsx/(Sh_lambda*f)
Sh_beta=2*np.pi/Sh_lambda
#結合に寄与するのは半波長コイルと仮定
l_cp=Sh_lambda/2
l_sh=0.9
l_sh=np.linspace(0.5, 3, num=25)
L1=Sh_L*l_cp
C1=Sh_C*l_cp
Cs=1/(L1*w**2)
#対象モデルでは半シート分のシート抵抗を利用
#定在波モデルでは1シート分のシート抵抗を利用
r1=1.0*l_sh*w*Sh_L/Q1
#カプラ側が2
L2=6.59*10**(-6)
# L2=np.linspace(0.1*10**(-6),1*10**(-6), num=20)
C2=1/(L2*w**2)
r2=w*L2/Q2
Lm=ke*np.sqrt(L1*L2)
Rl_opt=np.sqrt(r2**2+(r2*(w*Lm)**2)/(r1+Sh_Rsx))
# Rl_opt=np.linspace(1, 10, num=10)

# %%
def calcEtaAndZ_FromSource(w,l_FromSourceToCp,Sh_Rsx,Sh_beta,L1,Cs,r1,L2,C2,r2,Lm,RL):

    T=np.tan(Sh_beta*l_FromSourceToCp)
    Z1=r2+RL+1j*(w*L2-1/(w*C2))
    Z2=(w*Lm)**2/Z1
    Z3=r1+1j*w*L1+Z2+Sh_Rsx
    Z4=Sh_Rsx*(Z3+1j*Sh_Rsx*T)/(Sh_Rsx+1j*Z3*T)
    B=Z3-Z4
    Z5=Z4+B
    Zs=Z5-1j*(1/(w*Cs))

    Pr1=r1*(r2+RL)**2
    Pr2=r2*(w*Lm)**2
    PRsx=Sh_Rsx*(r2+RL)**2
    PRL=RL*(w*Lm)**2
    eta1=PRL/(Pr1+Pr2+PRsx+PRL)
    eta2=(PRL+PRsx)/(Pr1+Pr2+PRsx+PRL)

    # print("Z2 = " + str(Z2))
    # print("Z3 = " + str(Z3))
    # print("Z4 = " + str(Z4))
    # print("B = " + str(B))
    # print("Zs = " + str(Zs))
    # print("eta1 = " + str(eta1))
    # print("eta2 = " + str(eta2))

    return

# %%
def calcEtaAndZ_fromRL(w,Sh_Rsx,L1,r1,L2,r2,Lm):

    Z1=r1+2*Sh_Rsx+1j*w*L1
    Z2=((w*Lm)**2)/Z1
    Z3=r1+1j*w*L2+Z2
    Xc=Z3.imag
    Ccoil=1/(w*Xc)
    # Z4=Z3-1/(w*Ccoil)
    Z4=Z3.real

    Pr1=r1*((2*Sh_Rsx+r2)**2 + (w*Lm)**2)
    Pr2=r2*(w*Lm)**2
    Pz1=Sh_Rsx*(w*Lm)**2
    Pz2=Sh_Rsx*(w*Lm)**2
    eta=Pz1/(Pr1+Pr2+Pz1+Pz2)

    # print("Ccoil = " + str(Ccoil))
    # print("Z4_coil = " + str(Z4))
    # print("eta_fromCoil = " + str(eta))

    # with open('/mnt/c/Proj-2DWPT/calcForCoilSheet/eta_RL.csv', 'a') as f:
    #     writer = csv.writer(f)
    #     writer.writerow(eta)
    #     writer.writerow(Z4)
    #     writer.writerow(ke)

    return

# %%
def StandingWaveFromSource(w,Sh_Rsx,L1,r1,L2,r2,Lm):

    RL=(w*Lm)**2/(Sh_Rsx-r1)-r2
    RLopt_s=np.sqrt(r2**2 + r2*(w*Lm)**2/r1)
    Z_fromRL=(w*Lm)**2/(Sh_Rsx+r1)+r2

    Z4d=RL
    Z3d=Z4d+r2
    Z2d=(w*Lm)**2/Z3d
    Z1d=Z2d+r1

    eta_stand=RL*(w*Lm)**2/(r1*(RL*r2)**2+r2*(w*Lm)**2+RL*(w*Lm)**2)

    print("RL = " + str(RL))
    print("RLopt = " + str(RLopt_s))
    print("Z_fromRL = " + str(Z_fromRL))
    print("Z1d = " + str(Z1d))
    print("eta_stand = " + str(eta_stand))

    # with open('/mnt/c/Proj-2DWPT/calcForCoilSheet/CalcStandingWave.csv', 'a') as f:
    #     writer = csv.writer(f)
    #     writer.writerow(RL)
    #     writer.writerow(RLopt_s)
    #     writer.writerow(Z_fromRL)
    #     writer.writerow(Z1d)
    #     writer.writerow(eta_stand)

    return

# %%
calcEtaAndZ_FromSource(w,l_FromSourceToCp,Sh_Rsx,Sh_beta,L1,Cs,r1,L2,C2,r2,Lm,Rl_opt)
calcEtaAndZ_fromRL(w,Sh_Rsx,L1,r1,L2,r2,Lm)
StandingWaveFromSource(w,Sh_Rsx,L1,r1,L2,r2,Lm)
# %%
