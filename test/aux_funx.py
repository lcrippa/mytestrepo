import numpy as np
from edipy2 import global_env as ed
from mpi4py import MPI
import os,sys


def superconductive_zeta(warray,xmu,Sigma_all,axis):
    Ntot = np.shape(Sigma_all)[2]
    Lfreq = np.shape(warray)[0]
    zi = np.zeros((2,2,Ntot,Ntot,Lfreq),dtype=complex) #[Nambu,Nambu,Nso,Nso,Freq]
    
    if (axis == "m"):
        zi[0,0,:,:,:] = (warray + xmu) * np.eye(Ntot)[..., None]   -   Sigma_all[0]
        zi[0,1,:,:,:] =                                               -   Sigma_all[1]
        zi[1,0,:,:,:] =                                               -   Sigma_all[1]
        zi[1,1,:,:,:] = (warray - xmu) * np.eye(Ntot)[..., None]   +   np.conj(Sigma_all[0]) #Coleman 14.163
    elif (axis == "r"):
        warray_bar =   warray[::-1,...]
        Sigma_bar  =   Sigma_all[0,...,::-1]
        zi[0,0,:,:,:] = (warray + xmu) * np.eye(Ntot)[..., None]             - Sigma_all[0]
        zi[0,1,:,:,:] =                                                         - Sigma_all[1]
        zi[1,0,:,:,:] =                                                         - Sigma_all[1]
        zi[1,1,:,:,:] = -np.conj(warray_bar + xmu) * np.eye(Ntot)[..., None] + np.conj(Sigma_bar)
    return zi
    


def get_gloc(warray,xmu,Hk,Sigma_all,axis):
    '''
    Z has dimension  [Nambu,Nambu,Nso,Nso,Nfreq]
    Hk has dimension [Nnambu,Nk,Nso,Nso]
    Gk has dimension [Nk,Nnambu,Nso,Nso]
    Gmatrix has dimension [Nnambu*Ntot,Nnambu*Ntot,Nfreq]
    Returns an object of dimension [2,Ntot,Nfreq]
    '''

    try:
        import mpi4py
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        mpiflag = True
    except:
        mpiflag = False
        rank = 0
        size = 1

    master = (rank==0)
    
    if(master):
        print("Calculating local G axis "+axis+":")
    
    Z = superconductive_zeta(warray,xmu,Sigma_all,axis)
    Ntot = np.shape(Sigma_all)[2]
    Nfreq = np.shape(Z)[-1]
    Nk = np.shape(Hk)[1]
 
    
    if Nk >= Nfreq: 
        base = int(Nk // size)
        leftover = int(Nk % size)
        chunks = np.ones(size,dtype=int)* base
        chunks[:leftover] += 1
        offsets = np.zeros(size,dtype=int)
        offsets[1:]=np.cumsum(chunks)[:-1]
        ilow = offsets[rank]
        ihigh = ilow + chunks[rank]
        #print(rank,ilow,ihigh,chunks[rank])      
            
        Gtmp = np.zeros((chunks[rank],2*Ntot,2*Ntot,Nfreq),dtype=complex)
            
        Gtmp[:,0:Ntot,0:Ntot,:]             = Z[0,0,:,:,:] - Hk[0,ilow:ihigh,:,:,None]
        Gtmp[:,0:Ntot,Ntot:2*Ntot,:]        = Z[0,1,:,:,:] 
        Gtmp[:,Ntot:2*Ntot,0:Ntot,:]        = Z[1,0,:,:,:]
        Gtmp[:,Ntot:2*Ntot,Ntot:2*Ntot,:]   = Z[1,1,:,:,:] - Hk[1,ilow:ihigh,:,:,None]
       
        Gtmp = Gtmp.transpose(0, 3, 1, 2)
        Gtmp = np.linalg.inv(Gtmp)
        Gtmp = Gtmp.transpose(0, 2, 3, 1)
        if (mpiflag):
            Gtmp = np.ascontiguousarray(np.sum(Gtmp,axis=0)/Nk)
            Gloc = np.zeros_like(Gtmp)
            comm.Allreduce(Gtmp, Gloc, op=MPI.SUM)
        else:
            Gloc = np.sum(Gtmp,axis=0)/Nk
    else:
        base = int(Nfreq // size)
        leftover = int(Nfreq % size)
        chunks = np.ones(size,dtype=int)* base
        chunks[:leftover] += 1
        offsets = np.zeros(size,dtype=int)
        offsets[1:]=np.cumsum(chunks)[:-1]
        ilow = offsets[rank]
        ihigh = ilow + chunks[rank]
        #print(rank,ilow,ihigh,chunks[rank])      
        
        Gtmp = np.zeros((Nk,2*Ntot,2*Ntot,Nfreq),dtype=complex)
        Gloc = np.zeros((2*Ntot,2*Ntot,Nfreq),dtype=complex)
            
        Gtmp[:,0:Ntot,0:Ntot,ilow:ihigh]             = Z[0,0,:,:,ilow:ihigh] - Hk[0,:,:,:,None]
        Gtmp[:,0:Ntot,Ntot:2*Ntot,ilow:ihigh]        = Z[0,1,:,:,ilow:ihigh] 
        Gtmp[:,Ntot:2*Ntot,0:Ntot,ilow:ihigh]        = Z[1,0,:,:,ilow:ihigh]
        Gtmp[:,Ntot:2*Ntot,Ntot:2*Ntot,ilow:ihigh]   = Z[1,1,:,:,ilow:ihigh] - Hk[1,:,:,:,None]
       
        Gtmp = Gtmp.transpose(0, 3, 1, 2)
        Gtmp[:,ilow:ihigh,:,:] = np.linalg.inv(Gtmp[:,ilow:ihigh,:,:])
        Gtmp = Gtmp.transpose(0, 2, 3, 1)    
        
        if (mpiflag):
            Gtmp = np.ascontiguousarray(np.sum(Gtmp,axis=0)/Nk)    
            comm.Allreduce(Gtmp, Gloc, op=MPI.SUM)
        else:
            Gloc = np.sum(Gtmp,axis=0)/Nk

    return np.stack((Gloc[:Ntot,:Ntot,:], Gloc[:Ntot,Ntot:,:]), axis=0)

def dmft_weiss_field(G,Sigma):
    '''
    complex(8),dimension(2,Ntot,Ntot,Nfreq)              :: G
    complex(8),dimension(2,Ntot,Ntot,Nfreq)              :: Sigma
    '''
    try:
        import mpi4py
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        mpiflag = True
    except:
        mpiflag = False
        rank = 0
        size = 1

    master = (rank==0)
    
    if(master):
        print("Calculating Weiss Field")

    Ntot=np.shape(G)[1]
    Nfreq=np.shape(G)[-1]

    Gnambu = np.zeros((2*Ntot,2*Ntot,Nfreq),dtype=complex)
    Snambu = np.zeros((2*Ntot,2*Ntot,Nfreq),dtype=complex)

    Snambu[:Ntot,:Ntot,:]      = Sigma[0,:,:,:]
    Snambu[:Ntot,Ntot:,:]      = Sigma[1,:,:,:]
    Snambu[Ntot:,:Ntot,:]      = Sigma[1,:,:,:]
    Snambu[Ntot:,Ntot:,:]      = -np.conj(Sigma[0,:,:,:])
    
    Gnambu[:Ntot,:Ntot,:]      = G[0,:,:,:]
    Gnambu[:Ntot,Ntot:,:]      = G[1,:,:,:]
    Gnambu[Ntot:,:Ntot,:]      = G[1,:,:,:]
    Gnambu[Ntot:,Ntot:,:]      = -np.conj(G[0,:,:,:])

    Weiss = np.linalg.inv(np.linalg.inv(Gnambu.transpose(2, 0, 1)) + Snambu.transpose(2, 0, 1)).transpose(1, 2, 0)
    
    return np.stack((Weiss[:Ntot,:Ntot,:], Weiss[:Ntot,Ntot:,:]), axis=0)


def nn2so(mat):
    return mat.transpose(0, 2, 1, 3).reshape((ed.Nspin*ed.Norb, ed.Nspin*ed.Norb))

def so2nn(mat):
    return mat.reshape((ed.Nspin, ed.Norb, ed.Nspin, ed.Norb)).transpose(0, 2, 1, 3)
