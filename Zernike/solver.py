import numpy as np


def cg_solve(b,Amul,x0,epsilon=1e-6,nstepmax=500,Mdiv=None,args=None,verbose=True):

  x = x0.copy()
  nn=0
  r = b - Amul(x,*args)
  if (Mdiv is not None):
    z=Mdiv(r)
  else:
    z=r*1.
  p = z*1.
  rzold = np.vdot(z,r)
  while (nn<nstepmax):
    Ap = Amul(p,*args)
    alpha = rzold / np.vdot(p,Ap)
    x = x + alpha*p
    r = r - alpha*Ap
    print ("|r|=%f"%np.linalg.norm(r))

    if (cg_crit(r,epsilon)):
      print ("|r|=%f"%np.linalg.norm(r))
      return x
    if (Mdiv is not None):
      z=Mdiv(r)
    else:
      z=r*1.
    rz = np.vdot(r,z)
    beta = rz/rzold
    p = z + beta*p
    rzold = rz*1.
    nn = nn+1
  
def cg_crit(r,epsilon):

  normr = np.linalg.norm(r)
  if (normr < epsilon):
    return True
  else:
    return False
