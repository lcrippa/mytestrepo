import numpy as np
from edipy2 import global_env as ed
import mpi4py
from mpi4py import MPI
import os,sys
from aux_funx import *

#INIT MPI 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print("I am process",rank,"of",comm.Get_size())
master = (rank==0)

#Functions


def generate_kgrid(Nk):
    b1=2*np.pi*np.array([1.0,0.0])
    b2=2*np.pi*np.array([0.0,1.0])
    n1, n2 = np.meshgrid(np.arange(Nk), np.arange(Nk))
    n1=n1/Nk
    n2=n2/Nk
    gridout = np.stack([n1.ravel(), n2.ravel()], axis=-1)
    return np.dot(gridout,[b1,b2])

def h_square2d(k,t):
  return -2*t*(np.cos(k[...,0,np.newaxis,np.newaxis])+np.cos(k[...,1,np.newaxis,np.newaxis]))*np.eye(ed.Norb)

def test_dos(hk,plot=False):
    result=np.histogram(hk[:,0,0], bins=100, density=True) #normalized to one, I will have to multiply
    dx=result[1][2]-result[1][1]
    np.savetxt('dos.dat', np.transpose([np.delete(np.add(result[1],dx/2),np.size(result[1])-1), result[0]]))
    
    if plot:
        plt.xlabel('E')
        plt.ylabel("D(E)")
        plt.xlim(-3,3)
        plt.ylim(0,1)
        plt.plot(np.delete(np.add(result[1],dx/2),np.size(result[1])-1), result[0])
        plt.show()
    

#READ ED INPUT:
ed.read_input("inputAHM.conf")


#Parameters
Le      = 1000
wmixing = 0.5
wband   = 1.0
try:
    Nk = int(sys.argv[1])
except:
    Nk=20
    
print(ed.Uloc)
t_hop = 0.5

#BUILD frequency arrays and k grid:
wm = np.pi/ed.beta*(2*np.arange(ed.Lmats)+1)
wr = np.linspace(ed.wini,ed.wfin,ed.Lreal,dtype=complex)

kgrid = generate_kgrid(Nk)


#Generate hk and hloc
Hk   = h_square2d(kgrid,t_hop)
HkNambu   = np.array([h_square2d(kgrid,t_hop),-np.conj(h_square2d(-kgrid,t_hop))])
Hloc = np.sum(Hk,axis=0)/Nk**2
Hloc = Hloc.astype(complex)


#Generate dos and plot it
test_dos(Hk)

#SETUP SOLVER
ed.set_hloc(Hloc)
Nb=ed.get_bath_dimension()
bath = ed.init_solver()
bath_prev = np.copy(bath)


#DMFT CYCLE
converged=False;iloop=0
while (not converged and iloop<ed.Nloop ):
    iloop=iloop+1
    print("DMFT-loop:",iloop,"/",ed.Nloop)

    #Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
    ed.solve(bath)
    
    Smats = np.array([ed.get_sigma(axis="m",typ="n"),ed.get_sigma(axis="m",typ="a")])
    Sreal = np.array([ed.build_sigma(wr+1j*ed.eps,typ="n"),ed.build_sigma(wr+1j*ed.eps,typ="a")])
    

    Gmats = get_gloc(wm*1j       ,ed.xmu,HkNambu,Smats,axis="m")
    Greal = get_gloc(wr+1j*ed.eps,ed.xmu,HkNambu,Sreal,axis="r")  
    Weiss = dmft_weiss_field(Gmats,Smats)
          
    #Print       
    if(rank==0):
        for ispin in range(ed.Nspin):
            for iorb in range(ed.Norb):
                io=iorb+ed.Norb*ispin
                np.savetxt('Gloc_l'+str(iorb+1)+'_s'+str(ispin+1)+'_iw.dat', np.transpose([wm,Gmats[0,io,io,:].imag,Gmats[0,io,io,:].real]))
                np.savetxt('Floc_l'+str(iorb+1)+'_s'+str(ispin+1)+'_iw.dat', np.transpose([wm,Gmats[1,io,io,:].imag,Gmats[1,io,io,:].real]))
                np.savetxt('testSigma_l'+str(iorb+1)+'_s'+str(ispin+1)+'_iw.dat', np.transpose([wm,Smats[0,io,io,:].imag,Smats[0,io,io,:].real]))
                np.savetxt('testSelf_l'+str(iorb+1)+'_s'+str(ispin+1)+'_iw.dat', np.transpose([wm,Smats[1,io,io,:].imag,Smats[1,io,io,:].real]))
                np.savetxt('testSigma_l'+str(iorb+1)+'_s'+str(ispin+1)+'_realw.dat', np.transpose([np.real(wr),Sreal[0,io,io,:].imag,Sreal[0,io,io,:].real]))
                np.savetxt('testSelf_l'+str(iorb+1)+'_s'+str(ispin+1)+'_realw.dat', np.transpose([np.real(wr),Sreal[1,io,io,:].imag,Sreal[1,io,io,:].real]))
                np.savetxt('Gloc_l'+str(iorb+1)+'_s'+str(ispin+1)+'_realw.dat', np.transpose([np.real(wr),Greal[0,io,io,:].imag,Greal[0,io,io,:].real]))
                np.savetxt('Floc_l'+str(iorb+1)+'_s'+str(ispin+1)+'_realw.dat', np.transpose([np.real(wr),Greal[1,io,io,:].imag,Greal[1,io,io,:].real]))
                np.savetxt('Weiss_l'+str(iorb+1)+'_s'+str(ispin+1)+'_iw.dat', np.transpose([wm,Weiss[0,io,io,:].imag,Weiss[1,io,io,:].real]))

    #Fit
    ed.chi2_fitgf(Weiss[0],Weiss[1],bath)
        
    if(iloop>1):
        bath = wmixing*bath + (1.0-wmixing)*bath_prev
    bath_prev=np.copy(bath)
 
    err,converged=ed.check_convergence(Weiss,ed.dmft_error)

ed.finalize_solver()
print("Done...")

