from PlasmaState import *
from numpy import *

ps = PlasmaState("test", 1)
ps.read("iter.nc")
a= ps["nbion"]
print(a)

ps["nbion"][0] = 'H'
print(ps["nbion"])
ps.store("out.nc")

nrho = 20
rho =linspace(0.0,1.0,num=nrho)
ti = ps.interp1d("ti",0,rho)
print(ti)
nbeami = ps.interp1d("nbeami",0,rho,0)
print(nbeami)

