# make: gfortran -c -fPIC vr.F90; gfortran -shared vr.o -o libvr.so
# use: import vrpyintr
#      Qdp = vrpyintr.remap1(alg,dp1,dp2,Qdp)
# test: python vrpyintr.py

import ctypes

def load_vr():
    return ctypes.cdll.LoadLibrary("libvr.so")

def remap1(alg, dp1, dp2, Qdp):
    n = len(dp1)
    assert(len(dp2) == n)
    assert(len(Qdp) == n)
    assert(type(alg) == int)
    assert(alg >= -1 and alg <= 3)
    
    c_int = ctypes.c_int
    c_double = ctypes.c_double
    ndoubles = c_double * n
    def convert(arr): return ndoubles(*list(arr))

    Qdpc = convert(Qdp)
    load_vr().remap1c(c_int(n), c_int(alg), convert(dp1), convert(dp2), Qdpc)
    return list(Qdpc)

def test():
    import numpy as npy, math, matplotlib.pyplot as pl
    def cc(x): return 0.5*(x[:-1] + x[1:])
    p2 = npy.linspace(0, 1, 21)
    p2m = cc(p2)
    dp2 = npy.diff(p2)
    p1 = npy.sin(0.5*math.pi*p2)
    p1m = cc(p1)
    dp1 = npy.diff(p1)
    Qdp1 = npy.sin(2.2*math.pi*p1m + 0.3)*dp1
    Qdp2 = remap1(2, dp1, dp2, Qdp1)
    pl.plot(p1m, Qdp1/dp1, "r-", p2m, Qdp2/dp2, "k-")
    pl.show()

if __name__ == '__main__':
    test()
