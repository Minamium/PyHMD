"""Independent NumPy references for cell/face averages of oblique MHD waves."""
import numpy as np

FIELDS=('rho','mx','my','mz','energy','Bx','By','Bz','entropy_density','internal_energy')

def frame(config):
    g=config.grid;kx=2*np.pi/(g.x1_max-g.x1_min);ky=2*np.pi/(g.x2_max-g.x2_min)
    k=np.hypot(kx,ky)
    mode=config.physics.wave_mode
    speed={'fast':2.,'alfven':1.,'slow':.5,'entropy':1.}[mode] if config.problem=='linear_wave' else 1.
    return kx,ky,k,kx/k,ky/k,speed

def eigenvector(mode):
    if mode=='entropy':return (1.,0.,0.,0.,0.,0.,0.)
    if mode=='alfven':return (0.,0.,-1/3,np.sqrt(2)/1.5,0.,1/3,-np.sqrt(2)/1.5)
    c=2. if mode=='fast' else .5;den=c*c-1
    return (1/c,1.,-np.sqrt(2)/den,-.5/den,1/c,np.sqrt(2)*c/den,.5*c/den)

def point(config,x,y,time):
    kx,ky,k,nx,ny,speed=frame(config)
    phase=kx*(x-config.grid.x1_min)+ky*(y-config.grid.x2_min)-k*speed*time
    amplitude=config.physics.wave_amplitude;s=amplitude*np.sin(phase);c=amplitude*np.cos(phase)
    rho=np.ones_like(s);pressure=np.full_like(s,.1);bn=1.;bt=s;bz=c;vn=np.zeros_like(s);vt=-s;vz=-c
    if config.problem=='linear_wave':
        dr,dvn,dvt,dvz,dp,dbt,dbz=eigenvector(config.physics.wave_mode)
        rho=1+dr*s;pressure=.6+dp*s;vn=(1 if config.physics.wave_mode=='entropy' else 0)+dvn*s
        vt=dvt*s;vz=dvz*s;bt=np.sqrt(2)+dbt*s;bz=.5+dbz*s
    vx=nx*vn-ny*vt;vy=ny*vn+nx*vt;bx=nx*bn-ny*bt;by=ny*bn+nx*bt
    energy=pressure/(config.physics.gamma-1)+.5*rho*(vx*vx+vy*vy+vz*vz)+.5*(bx*bx+by*by+bz*bz)
    return np.stack((rho,rho*vx,rho*vy,rho*vz,energy,bx,by,bz,pressure*rho**(1-config.physics.gamma),pressure/(config.physics.gamma-1)),axis=-1)

def references(config,time):
    g=config.grid;dx=(g.x1_max-g.x1_min)/g.nx1;dy=(g.x2_max-g.x2_min)/g.nx2
    xf=np.linspace(g.x1_min,g.x1_max,g.nx1+1);yf=np.linspace(g.x2_min,g.x2_max,g.nx2+1)
    xc=.5*(xf[:-1]+xf[1:]);yc=.5*(yf[:-1]+yf[1:]);x,y=np.meshgrid(xc,yc)
    nodes,weights=np.polynomial.legendre.leggauss(8)
    volume=np.zeros((g.nx2,g.nx1,10))
    for i,a in enumerate(nodes):
        for j,b in enumerate(nodes):volume+=.25*weights[i]*weights[j]*point(config,x+.5*dx*a,y+.5*dy*b,time)
    kx,ky,k,nx,ny,speed=frame(config);coefficient=1.;bt=0.
    if config.problem=='linear_wave':coefficient=eigenvector(config.physics.wave_mode)[5];bt=np.sqrt(2)
    xx,yy=np.meshgrid(xf,yc);phase=kx*(xx-g.x1_min)+ky*(yy-g.x2_min)-k*speed*time
    bx=nx-ny*bt-ny*config.physics.wave_amplitude*coefficient*np.sin(phase)*np.sinc(ky*dy/(2*np.pi))
    xx,yy=np.meshgrid(xc,yf);phase=kx*(xx-g.x1_min)+ky*(yy-g.x2_min)-k*speed*time
    by=ny+nx*bt+nx*config.physics.wave_amplitude*coefficient*np.sin(phase)*np.sinc(kx*dx/(2*np.pi))
    discrete=volume.copy();discrete[:,:,5]=.5*(bx[:,:-1]+bx[:,1:]);discrete[:,:,6]=.5*(by[:-1,:]+by[1:,:])
    return volume,discrete,bx,by

def assemble(shards):
    first=shards[0];_,_,_,ny,_,_,_,nx=first.layout;ng=first.ng
    state=np.empty((ny,nx,first.nv));psi=np.empty((ny+1,nx+1))
    for shard in shards:
        _,_,local,ny,tx,ty,start,_=shard.layout
        state[:,start:start+local]=np.frombuffer(shard.buffers['U'],dtype='<f8').reshape(ty,tx,shard.nv)[ng:ng+ny,ng:ng+local]
        psi[:,start:start+local+1]=np.frombuffer(shard.buffers['Psi'],dtype='<f8').reshape(ty+1,tx+1)[ng:ng+ny+1,ng:ng+local+1]
    return state,psi

def errors(actual,expected):
    delta=np.abs(actual-expected)
    axes=tuple(range(delta.ndim-1))
    l1=np.mean(delta,axis=axes);linf=np.max(delta,axis=axes)
    return {'L1':dict(zip(FIELDS,map(float,l1))),'Linf':dict(zip(FIELDS,map(float,linf))),
            'rms_L1_physical':float(np.linalg.norm(l1[:8])),'max_Linf_physical':float(np.max(linf[:8]))}
