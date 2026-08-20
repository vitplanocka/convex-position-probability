import numpy as np, math
def cross(a,b): return a[...,0]*b[...,1]-a[...,1]*b[...,0]
def hull4_area(P):  # (m,4,2) -> convex hull area (works for convex-quad and interior-point cases)
    P0,P1,P2,P3=P[:,0],P[:,1],P[:,2],P[:,3]
    a1=0.5*np.abs(cross(P2-P0,P3-P1)); a2=0.5*np.abs(cross(P3-P0,P2-P1)); a3=0.5*np.abs(cross(P1-P0,P3-P2))
    t0=0.5*np.abs(cross(P2-P1,P3-P1)); t1=0.5*np.abs(cross(P2-P0,P3-P0))
    t2=0.5*np.abs(cross(P1-P0,P3-P0)); t3=0.5*np.abs(cross(P1-P0,P2-P0))
    return np.maximum.reduce([a1,a2,a3,t0,t1,t2,t3])
T=np.array([[[0,0],[1,0],[1,1],[0,1]],[[0,0],[1,0],[0,1],[0.1,0.1]],[[0,0],[2,0],[0,2],[0.5,0.5]]],float)
print("selftest (expect 1.0, 0.5, 2.0):",hull4_area(T))
def mc(body,N,seed=7,batch=5_000_000):
    rng=np.random.default_rng(seed); s2=0.0; n=0
    while n<N:
        m=min(batch,N-n)
        if body=='triangle':
            u=rng.random((m,4,2)); f=u.sum(2)>1; u[f]=1-u[f]; P=u; V=0.5
        elif body=='square': P=rng.random((m,4,2)); V=1.0
        elif body=='disk':
            r=np.sqrt(rng.random((m,4))); t=rng.random((m,4))*2*np.pi
            P=np.stack([r*np.cos(t),r*np.sin(t)],-1); V=math.pi
        A=hull4_area(P)/V; s2+=np.sum(A**2); n+=m
    return s2/N
tgt={'triangle':181/4500,'square':859/27000,'disk':(2400*math.pi**2+31031)/(19200*math.pi**4)}
for b in ['triangle','square','disk']:
    v=mc(b,60_000_000); e=tgt[b]; print(f"{b:9s} MC E[A4^2]={v:.8f}  target={e:.8f}  rel.diff={abs(v-e)/e:.2e}")
