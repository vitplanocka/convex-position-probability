"""Direct convex-position Monte Carlo of P_6(K) over a family of NON-regular convex bodies,
to stress-test the OPEN n>=6 extremal conjecture  P_6(triangle)=91/900 <= P_6(K) <= P_6(disk).
Uniform sampling in an arbitrary convex polygon via area-weighted fan triangulation.
"""
import numpy as np, math
from itertools import combinations
P6_TRI=91/900; P6_DISK=1-(146400*math.pi**2-473473)/(11520*math.pi**4)

def sample_poly(V,m,rng):
    V=np.asarray(V,float); n=len(V)
    tris=[(V[0],V[i],V[i+1]) for i in range(1,n-1)]
    ar=np.array([0.5*abs(np.cross(b-a,c-a)) for a,b,c in tris]); ar/=ar.sum()
    idx=rng.choice(len(tris),size=m,p=ar)
    u=rng.random(m); v=rng.random(m); f=u+v>1; u[f]=1-u[f]; v[f]=1-v[f]
    A=np.array([t[0] for t in tris])[idx]; B=np.array([t[1] for t in tris])[idx]; C=np.array([t[2] for t in tris])[idx]
    return A+ (B-A)*u[:,None] + (C-A)*v[:,None]

def in_cp(P):  # P:(m,6,2) convex position iff no point inside a triangle of 3 others
    m,n,_=P.shape; ok=np.ones(m,bool)
    def o(a,b,c): return (b[...,0]-a[...,0])*(c[...,1]-a[...,1])-(b[...,1]-a[...,1])*(c[...,0]-a[...,0])
    for i in range(n):
        others=[j for j in range(n) if j!=i]; pi=P[:,i]
        for a,b,c in combinations(others,3):
            A,B,C=P[:,a],P[:,b],P[:,c]; d1=o(A,B,pi); d2=o(B,C,pi); d3=o(C,A,pi)
            inside=((d1>=0)&(d2>=0)&(d3>=0))|((d1<=0)&(d2<=0)&(d3<=0)); ok&=~inside
    return ok

def mc_p6(V,N=100_000_000,seed=1,batch=2_000_000):
    rng=np.random.default_rng(seed); hit=0; n=0
    while n<N:
        b=min(batch,N-n); P=sample_poly(V,b*6,rng).reshape(b,6,2)
        hit+=int(in_cp(P).sum()); n+=b
    p=hit/N; return p, math.sqrt(p*(1-p)/N)

def reg(m,push=None,f=1.0):
    V=[[math.cos(2*math.pi*k/m),math.sin(2*math.pi*k/m)] for k in range(m)]
    if push is not None: V[push][0]*=f; V[push][1]*=f
    return [tuple(x) for x in V]
def halfdisk(n=64): return [(math.cos(math.pi*k/n),math.sin(math.pi*k/n)) for k in range(n+1)]

BODIES={
 "triangle":[(0,0),(1,0),(0,1)],
 "square":[(0,0),(1,0),(1,1),(0,1)],
 "reg-pentagon":reg(5),"reg-hexagon":reg(6),
 "trapezoid":[(0,0),(2,0),(1.3,1),(0.2,1)],
 "kite":[(0,0),(1,0.8),(0.3,2.0),(-0.6,0.8)],
 "pent-pushed":reg(5,0,1.5),
 "tri-pushed(quad)":[(0,0),(1,0),(0.2,0.2),(0,1)],   # non-convex? check
 "half-disk":halfdisk(64),
 "thin-tri-sliver":[(0,0),(1,0),(0.5,0.03)],          # affine to triangle -> must be 91/900
 "L-ish-hex":reg(6,0,1.6),
}
if __name__=="__main__":
    import sys
    print(f"window: triangle {P6_TRI:.6f}  disk {P6_DISK:.6f}\n")
    res=[]
    for name in (sys.argv[1:] or list(BODIES)):
        # skip clearly non-convex test bodies
        p,se=mc_p6(BODIES[name])
        z_lo=(p-P6_TRI)/se; z_hi=(P6_DISK-p)/se
        flag=""
        if p<P6_TRI-3*se: flag="  <-- BELOW TRIANGLE?!"
        if p>P6_DISK+3*se: flag="  <-- ABOVE DISK?!"
        print(f"{name:20s} P_6={p:.6f} +-{se:.1e}   above tri by {p-P6_TRI:+.5f} ({z_lo:+.0f}s), below disk by {P6_DISK-p:+.5f} ({z_hi:+.0f}s){flag}")
        res.append((name,p,se))
