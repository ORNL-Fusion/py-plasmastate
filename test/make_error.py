from PlasmaState import PlasmaState
import numpy as np
import sys
ps = PlasmaState("test", 1)

nrho = 51
nth = 101
ps["nrho_eq"] = nrho
ps["nth_eq"] = nth
ps["nrho_eq_geo"] = nrho
ps["nr"] = 0
ps["nz"] = 0

ps.alloc()

rho =np.linspace(0.0, 1.0, num=nrho)
ps["rho_eq"] = rho
ps["th_eq"] = np.linspace(0.0, 2.0*np.pi, num=nth)
ps["rho_eq_geo"] = rho

fn_geqdsk = sys.argv[1]

bdy_crat=1.e-6
kcur_option = 1
rho_curbrk = 0.9

ps["EQ_Code_Info"] = 'efit'
ps.updateEquilibrium(fn_geqdsk, bdy_crat, kcur_option, rho_curbrk)
ps.deriveMhdEq('everything')

ps.store("test.nc")

