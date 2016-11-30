import numpy as np
import math
from scipy import optimize as opt
from matplotlib import pyplot as plt


f = open('log.indent')                      #Opens the file
fw = open('Ni110.dat','w')
fw1 = open('Ni110More.dat','w')

varBool = False                            #Just a bool variable Load
firstCase = True                            #Another bool variable
secondCase = True                            #Another bool variable

h = []
P = []
t = []
frame = []
tempL = []
stepL = []
dt = 0.001 #ps
NStep = 1000
dtStep = 1.0
latticConst = 3.520
xlat = latticConst*np.sqrt(2)
ylat = latticConst*np.sqrt(2)
ev_AA = 1.602e-3 #ev/Å to uN
for line in f:
    if (line == 'Step Temp TotEng PotEng KinEng Press Pxx Pyy Pzz Lx Ly Lz v_y f_4 f_4[1] f_4[2] f_4[3] \n'):            #Checks line
        varBool = True
        continue
    if varBool:
        if (line.split()[0] == "Loop"):
            varBool = False
            if int(line.split()[8]) == 20000: #20.000 is chosen hold time
                continue
            elif firstCase:
                firstCase = False
                continue
            else:
                break
        step = int(line.split()[0])        #Typecast step# to int
        temp = float(line.split()[1])        #Typecast temp to float
        toteng = float(line.split()[2])       #Typecast total Energy to float
        poteng = float(line.split()[3])       #Typecast Pressure to float
        kineng = float(line.split()[4])       #Typecast volume to float
        Press = float(line.split()[5])       #Typecast pressure in y to float
        Pxx = float(line.split()[6])       #Typecast box size to float
        Pyy = float(line.split()[7])       #Typecast box size to float
        Pzz = float(line.split()[8])       #Typecast box size to float
        Lx = float(line.split()[9])       #Typecast box size to float
        Ly = float(line.split()[10])       #Typecast box size to float
        Lz = float(line.split()[11])       #Typecast box size to float
        v_y = float(line.split()[12])       #Typecast box size to float
        f_4 = line.split()[13:17]       #Typecast box size to float
        time = step*0.001
        if time == 0:
            continue
        frame.append(step/NStep)
        h.append(v_y)
        P.append(float(f_4[2])*ev_AA) #uN
        t.append(time)
        tempL.append(temp)
        stepL.append(step)

h = [-(i-max(h))*1e-1 for i in h] #nm
fw.write("h"+" "+"p"+" "+"t"+"\n")
for i,j,k in zip(h,P,t):
    fw.write(str(i)+" "+str(j)+" "+str(k)+"\n")
fw1.write("h"+" "+"p"+" "+"t"+" "+"temp"+" "+"step"+"\n")
for i,j,k,l,m,n in zip(h,P,t,tempL,stepL,frame):
    fw1.write(str(i)+" "+str(j)+" "+str(k)+" "+str(l)+" "+str(m)+" "+str(n)+"\n")
exit()
#plt.plot(h,P)
plt.plot(stepL, [i/max(P) for i in P])
plt.plot(stepL, [i/max(tempL) for i in tempL])
plt.show()
exit()
