"""Two-chord Blaschke-Petkantschin route to E[A_4^2 & convex] for the DISK (validation).
E[A_4^2 & convex] = (3/4)(1/V^4) int_{l1,l2} sin^2(psi) G(chord1;sigma1) G(chord2;sigma2) 1{X in K} dl1 dl2,
G(-h,h,sigma) = ((2h)^5 - (h-sigma)^5 - (sigma+h)^5)/10   [ = int int |s-s'|^3 1{straddle sigma} ].
Disk: fix phi1=0 (factor pi by rotation), integrate psi in (0,pi), p1,p2 in (-1,1); X = intersection.
Target (Marckert): (2400 pi^2 + 19019)/(19200 pi^4) = 0.0228343519...
"""
import numpy as np
PI=np.pi
def G(h,sig):
    return ((2*h)**5-(h-sig)**5-(sig+h)**5)/10
def integrand_grid(Npsi=700,Np=700):
    psi=(np.arange(Npsi)+0.5)/Npsi*PI; dpsi=PI/Npsi
    p=(np.arange(Np)+0.5)/Np*2-1; dp=2.0/Np
    tot=0.0
    s2=np.sin(psi)**2; sn=np.sin(psi); cs=np.cos(psi)
    P1,P2=np.meshgrid(p,p,indexing='ij')       # p1,p2 grid
    h1=np.sqrt(np.clip(1-P1*P1,0,None)); h2=np.sqrt(np.clip(1-P2*P2,0,None))
    for k in range(Npsi):
        yX=(P2-P1*cs[k])/sn[k]
        inside=(P1*P1+yX*yX)<=1.0
        sig1=yX
        sig2=-P1*sn[k]+yX*cs[k]
        val=s2[k]*G(h1,sig1)*G(h2,sig2)
        tot+=np.sum(np.where(inside,val,0.0))*dp*dp*dpsi
    return tot
I=integrand_grid()
EA=0.75*(1.0/PI**6)*PI*I
target=(2400*PI**2+19019)/(19200*PI**4)
print(f"two-chord grid E[A4^2&conv](disk) = {EA:.8f}")
print(f"target (Marckert)                 = {target:.8f}")
print(f"rel.diff = {abs(EA-target)/target:.2e}")
