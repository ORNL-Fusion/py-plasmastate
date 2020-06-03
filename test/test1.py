#!/usr/bin/env python
# coding: utf-8

from PlasmaState import *
from numpy import *

ps = PlasmaState("test", 1)

nspec_rfmin = 1
nspec_th = 2
nrho = 51
nth = 101

ps["nspec_th"] = nspec_th 
ps["nspec_rfmin"] = nspec_rfmin
ps["nrho"] =  nrho
ps["nrho_eq"] = nrho
ps["nth_eq"] = nth
ps["nrho_eq_geo"] = nrho
ps["nr"] = 0
ps["nz"] = 0

version = ps.getVersionId()
label   = "plasma_state_test"

ps["global_label"] = label + " " + version
ps["runid"] = "Test"
ps["tokamak_id"] = "Dummy"
ps["shot_number"] = 12345

time0 = 0.01
ps["t0"] = time0
ps["t1"] = time0

ps.alloc()

print('allocated')

ps.setThermalSpecies(-1, -1, 0)
ps.setThermalSpecies(1, 1, 2)
izImp = 6
ps.setThermalSpecies(izImp, izImp, 0)

ps.setRFMinoritySpecies(1, 1, 1)
ps.finishSpecies()

rho =linspace(0.0,1.0,num=nrho)
ps["rho"] = rho
ps["rho_eq"] = rho
ps["th_eq"] = linspace(0.0,2.0*pi,num=nth)
ps["rho_eq_geo"] = rho

print(ps["rho"])

yn = 1.0-rho**2
yc = [ 0.5*(yn[i]+yn[i+1]) for i in range(nrho-1) ]
ps["ti"] = yc

ps.store("out.nc");

nrho = 20
rho =linspace(0.0,1.0,num=nrho)
tmp = ps.interp1d("ti",0,rho)
print(tmp)
