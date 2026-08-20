import numpy as np, math
def cross(a,b): return a[...,0]*b[...,1]-a[...,1]*b[...,0]
def seg_cross(P1,P3,P2,P4):  # do segments P1P3 and P2P4 intersect?
    d1=np.sign(cross(P3-P1,P2-P1)); d2=np.sign(cross(P3-P1,P4-P1))
    d3=np.sign(cross(P4-P2,P1-P2)); d4=np.sign(cross(P4-P2,P3-P2))
    return (d1*d2<0)&(d3*d4<0)
def hull4_area(P):
    P0,P1,P2,P3=P[:,0],P[:,1],P[:,2],P[:,3]
    a1=0.5*np.abs(cross(P2-P0,P3-P1)); a2=0.5*np.abs(cross(P3-P0,P2-P1)); a3=0.5*np.abs(cross(P1-P0,P3-P2))
    t0=0.5*np.abs(cross(P2-P1,P3-P1)); t1=0.5*np.abs(cross(P2-P0,P3-P0))
    t2=0.5*np.abs(cross(P1-P0,P3-P0)); t3=0.5*np.abs(cross(P1-P0,P2-P0))
    return np.maximum.reduce([a1,a2,a3,t0,t1,t2,t3])
def convex4(P):  # all 4 points are hull vertices <=> exactly one split crosses
    P0,P1,P2,P3=P[:,0],P[:,1],P[:,2],P[:,3]
    c1=seg_cross(P0,P2,P1,P3); c2=seg_cross(P0,P1,P2,P3); c3=seg_cross(P0,P3,P1,P2)
    return c1|c2|c3
def mc(body,N=60_000_000,seed=11,batch=5_000_000):
    rng=np.random.default_rng(seed); sA=0.0; sD=0.0; n=0
    while n<N:
        m=min(batch,N-n)
        if body=='triangle':
            u=rng.random((m,4,2)); f=u.sum(2)>1; u[f]=1-u[f]; P=u; V=0.5
        elif body=='square': P=rng.random((m,4,2)); V=1.0
        else:
            r=np.sqrt(rng.random((m,4))); t=rng.random((m,4))*2*np.pi
            P=np.stack([r*np.cos(t),r*np.sin(t)],-1); V=math.pi
        conv=convex4(P)
        A=hull4_area(P)/V
        sA+=np.sum((A**2)*conv)                       # E[A4^2 & convex]
        D=cross(P[:,2]-P[:,0],P[:,3]-P[:,1])/V        # det(P3-P1,P4-P2)/V
        cr=seg_cross(P[:,0],P[:,2],P[:,1],P[:,3])
        sD+=np.sum((D**2)*cr)                          # E[D^2 1{P1P3 x P2P4}]
        n+=m
    return sA/N, sD/N
tgt={'triangle':119/4500,'square':1307/54000,'disk':(2400*math.pi**2+19019)/(19200*math.pi**4)}
for b in ['triangle','square','disk']:
    EAconv, ED2 = mc(b); pred=0.75*ED2; e=tgt[b]
    print(f"{b:9s} E[A4^2&conv](direct)={EAconv:.7f}  (3/4)E[D^2 1cross]={pred:.7f}  exact={e:.7f}")
