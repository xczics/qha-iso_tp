import ctypes as ct
import numpy as np
#gfortran -shared Fdistance.f90 -o Fdistance.so 
import os
fpath=os.path.dirname(__file__)
#trycompile=os.system("gfortran -shared beta.f90 -o Fbeta.so")
if not os.path.isfile(fpath+'/Fbeta.so'):
    trycompile=os.system("gfortran -shared %s/beta.f90 -o %s/Fbeta.so"%(fpath,fpath))
    if trycompile:
        print("无动态链接文件，且尝试编译失败，程序退出！请检查gfortran是否可用，或自行编译beta.f90 为 Fbeta.so并粘贴到目录flib下")
        exit()
fortlib = ct.CDLL(fpath+'/Fbeta.so')
_GetLnbeta=fortlib.GetLnbeta
#GetLnbeta (Hfreq,Lfreq,Weights,Nq,Nfreq,An,T,NT,Lnbeta)
#integer(c_int), intent(in) :: Nq, Nfreq,An,NT
#real(c_double), intent(in) :: Lfreq(Nq,Nfreq), Hfreq(Nq,Nfreq), Weights(Nq)
#real(c_double), intent(in) :: T(NT)
#real(c_double), intent(inout) :: Lnbeta(NT)

_GetLnbeta.argtypes=[ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.c_int,ct.c_int,ct.c_int,ct.POINTER(ct.c_double),ct.c_int,ct.POINTER(ct.c_double)]

def GetLnbeta(Hfreq,Lfreq,Weights,An,T):
    T=np.array(T,dtype=ct.c_double,order='F')
    Hfreq=np.array(Hfreq,dtype=ct.c_double,order='F')
    Lfreq=np.array(Lfreq,dtype=ct.c_double,order='F')
    Weights=np.array(Weights,dtype=ct.c_double,order='F')
    Nq=Hfreq.shape[0]
    Nfreq=Hfreq.shape[1]
    NT=len(T)
    Lnbeta=np.zeros(NT)
    #print(Hfreq)
    _GetLnbeta(Hfreq.ctypes.data_as(ct.POINTER(ct.c_double)),
               Lfreq.ctypes.data_as(ct.POINTER(ct.c_double)),
               Weights.ctypes.data_as(ct.POINTER(ct.c_double)),
               ct.c_int(Nq),ct.c_int(Nfreq),ct.c_int(An),
               T.ctypes.data_as(ct.POINTER(ct.c_double)),
               ct.c_int(NT),
               Lnbeta.ctypes.data_as(ct.POINTER(ct.c_double)))
    return Lnbeta
