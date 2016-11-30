import numpy as np
import math
from scipy import optimize as opt
from matplotlib import pyplot as plt
from numpy import inf as inf


#Cu = open('Cu.dat')                      #Opens the file
#Cu111 = open('Cu111More.dat')                      #Opens the file
#Ni010 = open('Ni010More.dat')                      #Opens the file
Ni110 = open('Ni110More.dat')                      #Opens the file
#Fq = open('FQ.dat')


def pressureLoad(h,alpha,n):
    return alpha*h**n
def pressureUnLoad(h,alpha,m,h_f):
    return alpha*(h-h_f)**m
def stiffnessUnLoad(h, alpha, m, h_f):
    return m*alpha*(h-h_f)**(m-1)
def geth_c(h_max, eps, P_max, S):
    return h_max-eps*P_max/S

def findIndex(Pylist, threshold):
    return [i for i in range(len(Pylist)) if Pylist[i] > threshold*max(Pylist)][0]

def findFirst(Pylist, cond):
    for i in Pylist:
        print(i)
    return 0

def Berkovich(h_c):
    return 24.56*h_c**2-5.84e4*h_c+5.73e6*h_c**(0.5)-6.92e7*h_c**(0.25)+1.9137e8*h_c**(1/8)-1.30e8*h_c**(1/16)
def Berkovich_perfect(h_c):
    return 24.56*h_c**2#-5.84e4*h_c+5.715e6*h_c**(0.5)-6.916e7*h_c**(0.25)+1.9137e8*h_c**(1/8)-(1.30e8)*(h_c**(1/16))
def Cylinder(h_c, uni):
    if uni == 'AA':
        R = 5
        Z = 4
    elif uni == 'nm':
        R = 0.5
        Z = 0.4
    elif uni == 'um':
        R = 5*1e-4
        Z = 4*1e-4
    else:
        print('No valid units given, either AA, nm or um')
        return False
    return 2*Z*np.sqrt(h_c*(2*R-h_c)) # Units uni

#exit()

def FitCurve(h, P, fun, BL = (-inf, inf)): #FC = fittedCurve, Cov = Covariance
    if fun.__code__.co_argcount == 4:
        FC, cov=opt.curve_fit(fun, h, P, bounds = BL)
    else:
        FC, cov=opt.curve_fit(fun, h, P, bounds = BL)
    return FC, cov

def plotAFun():
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

def routine1(f, Llim, ULlim):
    header = next(f)
    lines = f.readlines()
    depth = [float(i.split()[0]) for i in lines] #[h] = nm
    load = [float(i.split()[1]) for i in lines] #[P] = uN
    time = [float(i.split()[2]) for i in lines] #[t] = s
    temp = [float(i.split()[3]) for i in lines] #[t] = s
    step = [float(i.split()[4]) for i in lines] #[t] = s
    frame = [float(i.split()[5]) for i in lines] #[t] = s
    
    #StartIndex, EndIndex
    SILoad = findIndex(load, Llim[0])
    EIUnLoad = len(load)-findIndex(load[::-1], ULlim[0])
    EILoad = findIndex(load, Llim[1])
    SIUnLoad = len(load)-findIndex(load[::-1], ULlim[1])
    
    #Curve Fit
    (alphaL, n), Lcovlist = FitCurve(depth[SILoad:EILoad],load[SILoad:EILoad], pressureLoad)
    BUL = ([0,0,0.0],[inf,6,1.3])
    (alphaUL, m, h_f), ULcovlist = FitCurve(depth[SIUnLoad:EIUnLoad],load[SIUnLoad:EIUnLoad], pressureUnLoad, BUL)

    #Plotting
    plt.figure(figsize = (14,10), dpi=160)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.rc('axes', labelsize=18)
    plt.rc('legend', fontsize=18)

    
    plt.plot(depth, load,'--')
    plt.plot(depth[SILoad:EILoad], load[SILoad:EILoad], label='load part')
    plt.plot(depth[SIUnLoad:EIUnLoad], load[SIUnLoad:EIUnLoad], label='unload part')
 #   plt.show(),exit()
 #   plt.show()
#    (alphaL, n), Lcovlist = FitCurve(depth[SILoad:EILoad],load[SILoad:EILoad], pressureLoad)
#    (alphaUL, m, h_f), ULcovlist = FitCurve(depth[SIUnLoad:EIUnLoad],load[SIUnLoad:EIUnLoad], pressureUnLoad)
#    print('AlphaL = {0}\nAlphaUL = {1}\n n,m = {2},{3}\n h_f = {4}'.format(alphaL, alphaUL, n,m,h_f))
#    exit()
    plt.plot(depth[SILoad:EILoad], [pressureLoad(i, alphaL, n) for i in depth[SILoad:EILoad]], \
    label=r'$F = {:.2f}*h^{{{:.2f}}}$'.format(alphaL,n))
    plt.plot(depth[SIUnLoad:EIUnLoad], [pressureUnLoad(i, alphaUL, m, h_f) for i in depth[SIUnLoad:EIUnLoad]], \
    label=r'$F = {:.2f}*(h-{:.2f})^{{{:.2f}}}$'.format(alphaUL,h_f,m))
#    plt.show(),exit()
    
    plt.xlabel(r'Displacement $h$ [nm]')
    plt.ylabel(r'Load $F$ [$\mu$ N]')
    plt.title(r'Nano indentation on Nickel in [110]', fontsize = 26, y=1.04)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc='best')
    plt.show(),exit()
    plt.savefig('NanoI{}.png'.format(f.name.split(".")[0]), transparent=True)
    plt.close()


    exit()

    print('AlphaL = {0}\nAlphaUL = {1}\n n,m = {2},{3}\n h_f = {4}'.format(alphaL, alphaUL, n,m,h_f))
    P_max = max(load)
    h_max = max(depth)
    S = stiffnessUnLoad(h_max, alphaUL, m, h_f)
    h_c = geth_c(h_max, 0.75, P_max, S)
    A = Cylinder(h_c*1e-3,'um')
    E_r = np.pi**0.5/2/A**0.5*S
    H = P_max/A*1e-3
    print('\n\nThis is for {}'.format(header))
    print('S = \t{:.3f} GPa'.format(S))
    print('h_c = \t{:.3f} nm'.format(h_c))
    print('A = \t{:.3f} (uN)²'.format(A))
    print('E_r = \t{:.3f} GPa'.format(E_r))
    print('H = \t{:.3f} GPa'.format(H))
    return h_c



def main():
   # routine1(Ni010, [0.01, 0.6], [0.04, 0.5])
    routine1(Ni110, [0.01, 0.9], [0.04, 0.5])
    return 0
    

main()
