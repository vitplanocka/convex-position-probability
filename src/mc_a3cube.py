import numpy as np, math
def mc(body,N=200_000_000,seed=3,batch=5_000_000):
    rng=np.random.default_rng(seed); s3=0.0; n=0
    while n<N:
        m=min(batch,N-n)
        if body=='triangle':
            u=rng.random((m,3,2)); f=u.sum(2)>1; u[f]=1-u[f]; P=u; V=0.5
        elif body=='square':
            P=rng.random((m,3,2)); V=1.0
        elif body=='disk':
            r=np.sqrt(rng.random((m,3))); t=rng.random((m,3))*2*np.pi
            P=np.stack([r*np.cos(t),r*np.sin(t)],-1); V=math.pi
        d=(P[:,1,0]-P[:,0,0])*(P[:,2,1]-P[:,0,1])-(P[:,1,1]-P[:,0,1])*(P[:,2,0]-P[:,0,0])
        A=np.abs(d)/2/V
        s3+=np.sum(A**3); n+=m
    return s3/N
for b,exact in [('triangle',31/9000),('square',137/72000),('disk',1001/(6400*math.pi**4))]:
    v=mc(b); print(f"{b:9s} MC E[A3^3]={v:.8f}  exact={exact:.8f}  rel.diff={abs(v-exact)/exact:.2e}")
