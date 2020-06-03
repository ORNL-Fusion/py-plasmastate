from PlasmaState import *
ps = PlasmaState("ips",1)
ps.read("ips-state.nc")
print(ps["nbion"])
print(ps["nbi_src_name"])

