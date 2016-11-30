import numpy as np
import math
from scipy import optimize as opt
from matplotlib import pyplot as plt


#Open data for reading

#Cu = open('Cu.dat')                      #Opens the file
#Cu110 = open('Cu110.dat')                      #Opens the file
#Ni = open('Ni.dat')                      #Opens the file
#Ni110 = open('Ni110.dat')                      #Opens the file
Fq = open('FQ.dat')
Al = open('Al.dat')


#Create/overwrite file for writing
#Fw = open('X.dat','w')

def pressureLoad(h,alpha,n): #General formula for the loading relationship between load and displacement/depth
    return alpha*h**n
def pressureUnLoad(h,alpha,m,h_f):#General formula for the loading relationship between un-loading and displacement/depth
    return alpha*(h-h_f)**m
def stiffnessUnLoad(h, alpha, m, h_f): #Formula for the stiffness
    return m*alpha*(h-h_f)**(m-1)
def pressureHertz(h, R_h):  #Load using the Hertz formula. Only in the elastic region
    return R_h**0.5*h**(3/2)
def geth_c(h_max, eps, P_max, S):   #calculates h_c
    return h_max-eps*P_max/S

def findIndex(Pylist, threshold): #Used to find start and end of loading/unloading
    return [i for i in range(len(Pylist)) if Pylist[i] > threshold*max(Pylist)][0]

#Area function
def Berkovich(h_c): #h_c > 200
    return 24.56*h_c**2-5.84e4*h_c+5.73e6*h_c**(0.5)-6.92e7*h_c**(0.25)+1.9137e8*h_c**(1/8)-1.30e8*h_c**(1/16)
def Berkovich_perfect(h_c):
    return 24.56*h_c**2 # For h_c < 150

#exit()

def FitCurve(h, P, fun, maxfev_T = 800, endval=0): #FC = fittedParameters, Cov = Covariance
    FC, cov=opt.curve_fit(fun, h, P, maxfev = maxfev_T)
    return FC, cov

def plotAFun(): #Plots the areafunction to look at its validity
    lio = [30 +i/10 for i in range(int(200e1))]
    lia = [Berkovich(i) for i in lio]
    index = findIndex(lia,0)
    plt.figure(figsize = (14,10), dpi=160)
    plt.plot(lio, lia)
    plt.plot(lio[index], [0], 'x', label='h_c = {}'.format(lio[index]))
    plt.xlabel('h_c [nm]')
    plt.ylabel(r'A $\left[\frac{1}{nm^2}\right]$')
    plt.title('Berkovich Area Function')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.savefig('BerkovichAF.png', transparent=True)
    plt.close()


def routine1(f):
    header = next(f)
    lines = f.readlines()
    depth = [float(i.split()[0]) for i in lines] #[h] = nm
    load = [float(i.split()[1]) for i in lines] #[P] = uN
    time = [float(i.split()[2]) for i in lines] #[t] = s
#StartIndex, EndIndex
    SILoad = findIndex(load, 0.05)
    EIUnLoad = len(load)-findIndex(load[::-1], 0.05)
    EILoad = findIndex(load, 0.95)
    SIUnLoad = len(load)-findIndex(load[::-1], 0.95)

    (alphaL, n), Lcovlist = FitCurve(depth[SILoad:EILoad],load[SILoad:EILoad], pressureLoad)
    (alphaUL, m, h_f), ULcovlist = FitCurve(depth[SIUnLoad:EIUnLoad],load[SIUnLoad:EIUnLoad], pressureUnLoad)

    print('AlphaL = {0}\nAlphaUL = {1}\n n,m = {2},{3}\n h_f = {4}'.format(alphaL, alphaUL, n,m,h_f))
    P_max = max(load)
    h_max = max(depth)
    S = stiffnessUnLoad(h_max, alphaUL, m, h_f)
    h_c = geth_c(h_max, 0.75, P_max, S)
 #   A = (Berkovich(h_c) if h_c > 200 else Berkovich_perfect(h_c))*1e-6 #[A]=um²
 #   A = (Berkovich(h_c*1e1) if h_c > 200 else Berkovich_perfect(h_c*1e1))*1e-8 #[A]=um²
    A = Berkovich_perfect(h_c*1e-3)
    E_r = np.pi**0.5/2/A**0.5*S
    H = P_max/A*1e-3
    h_elast = h_max-h_f
    print('\n\nThis is for {}'.format(header))
    print('S = \t{:.3f} GPa'.format(S))
    print('h_c = \t{:.3f} nm'.format(h_c))
    print('h_e = \t{:.3f} nm'.format(h_elast))
    print('A = \t{:.3f} (uN)²'.format(A))
    print('E_r = \t{:.3f} GPa'.format(E_r))
    print('H = \t{:.3f} GPa'.format(H))
 #   h_elast = h_max-h_c

    EIHert = findIndex(depth, h_elast*0.9/max(depth))
    SIHert = [i for i in range(len(depth)) if depth[i] == 0][0]

    R, Hcov = FitCurve(depth[SIHert:EIHert], load[SIHert:EIHert], pressureHertz)
    R = R[0]
    R_val = R/(4/3*E_r)**2*1e6
    print('R = \t{:.3f} nm'.format(R_val))

    plt.figure(figsize = (14,10), dpi=160)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.rc('axes', labelsize=18)
    plt.rc('legend', fontsize=18)

    labUL = 'Load: $P = {:.2f}*(h-{:.2f})^{{{:.2f}}}$'.format(alphaUL, h_f,m)
    labL = r'UnLoad: $P = {:.2f}*h^{{{:.2f}}}$'.format(alphaL, n)
    labH = r'Hertz: $P = {:.2f}^{{\frac{{1}}{{2}}}}*h^{{\frac{{3}}{{2}}}}$'.format(R)
    lab = 'fused quartz'
    if f.name.split(".")[0] == 'Al':
        lab = 'aluminium'
    plt.plot(depth, load,'--')
    plt.plot(depth[SILoad:EILoad], load[SILoad:EILoad],'b', label='Experimental data')
    plt.plot(depth[SIUnLoad:EIUnLoad], load[SIUnLoad:EIUnLoad],'b')
    plt.plot(depth[SIHert:EIHert], load[SIHert:EIHert],'b')
    plt.plot(depth[SILoad:EILoad], [pressureLoad(i, alphaL, n) for i in depth[SILoad:EILoad]], label = labL)
    plt.plot(depth[SIUnLoad:EIUnLoad], [pressureUnLoad(i, alphaUL, m, h_f) for i in depth[SIUnLoad:EIUnLoad]], label=labUL)
    plt.plot(depth[SIHert:EIHert], [pressureHertz(i, R) for i in depth[SIHert:EIHert]], label=labH)
    plt.xlabel(r'Displacement $h$ [nm]')
    plt.ylabel(r'Load $P$ [$\mu N$]')
    plt.title(r'Nano indentation on {}'.format(lab), fontsize = 26, y=1.04)
    plt.ylim(min(load), max(load)+100)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc="best")
 #   plt.show()#, exit()
    return 0
    plt.savefig('NanoI{}.png'.format(f.name.split(".")[0]), transparent=True)
    plt.close()
    
    return 0




def main():
    routine1(Al)
    routine1(Fq)
#    routine1(Cu)
#    routine1(Cu110)
#    routine1(Ni)
#    routine1(Ni110)
#    plotAFun()
#    plt.show()

    return 0
    

main()
