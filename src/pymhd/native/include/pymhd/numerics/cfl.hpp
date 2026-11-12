#pragma once
#include "pymhd/numerics/thermal.hpp"

// The full record is copied between device/host and broadcast as bytes. Keep
// every byte initialized: an implicit tail pad triggers CUDA initcheck even
// though it is not used in the numerical reduction.
struct Rate {double adv=0,diff=0,bad=0,total=0;int cell=-1,reserved=0;};
static_assert(sizeof(double)==8&&sizeof(int)==4&&sizeof(Rate)==40,"CFL transfer record must have no implicit padding");
HD inline bool better_rate(double total,int cell,double previous_total,int previous_cell) {
    return total>previous_total||(total==previous_total&&cell>=0&&(previous_cell<0||cell<previous_cell));
}
HD inline Rate combine_rates(Rate left,Rate right) {
    double bad=fmax(left.bad,right.bad);
    if(better_rate(right.total,right.cell,left.total,left.cell))left=right;
    left.bad=bad;return left;
}
HD inline Rate rate_at(int k,const State* u,const double* a,MeshView m,Parameters p) {
    int i=k%m.nx+NG,j=k/m.nx+NG;State q=u[m.cell(i,j)];Rate r;
    r.cell=(j-NG)*m.global_nx+m.start+i-NG;
    for(int c=0;c<NV;++c)if(!finite_value(q[c]))r.bad=1;
    const double pressure=raw_pressure(q,m.radius(i,j,p),p);
    if(!(pressure>0) || !finite_value(pressure))r.bad=1;
    if(r.bad)return r;
    Primitive v=get_primitive(u,i,j,m,p);
    double dx=m.xf[i+1]-m.xf[i],dy=(m.yf[j+1]-m.yf[j])*(p.spherical?m.xc[i]:1.);
    Primitive wave=v;
    if(p.cfl_staggered) {
        // Opposing face fields may cancel in cell B, but are still used as
        // Riemann normal fields and by CT. Do not hide their wave speeds.
        double bl=field1(a,i,j,m,p),br=field1(a,i+1,j,m,p),bd=field2(a,i,j,m,p),bu=field2(a,i,j+1,m,p);
        if(!finite_value(bl)||!finite_value(br)||!finite_value(bd)||!finite_value(bu)){r.bad=1;return r;}
        wave[B1]=fmax(fabs(v[B1]),fmax(fabs(bl),fabs(br)));
        wave[B2]=fmax(fabs(v[B2]),fmax(fabs(bd),fabs(bu)));
    }
    r.adv=(fabs(v[M1])+fast_speed(wave,0,p))/dx+(fabs(v[M2])+fast_speed(wave,1,p))/dy;
    double ep=eta_cell_value(u,i,j,m,p)*fmax(1.,p.eta_model==0?1.:p.chi);
    r.diff=2*ep*(1/(dx*dx)+1/(dy*dy));r.total=r.adv+r.diff;return r;
}
#ifdef PYMHD_CUDA
__global__ void rate_reduce(const State* u,const double* a,Rate* results,MeshView m,Parameters p) {
    __shared__ double adv[128],diff[128],bad[128],total[128];__shared__ int cell[128];int k=blockIdx.x*128+threadIdx.x,t=threadIdx.x;
    Rate r;if(k<m.nx*m.ny)r=rate_at(k,u,a,m,p);adv[t]=r.adv;diff[t]=r.diff;bad[t]=r.bad;total[t]=r.total;cell[t]=r.cell;__syncthreads();
    for(int stride=64;stride;stride/=2){if(t<stride){
        if(better_rate(total[t+stride],cell[t+stride],total[t],cell[t])){adv[t]=adv[t+stride];diff[t]=diff[t+stride];total[t]=total[t+stride];cell[t]=cell[t+stride];}
        bad[t]=fmax(bad[t],bad[t+stride]);}__syncthreads();}
    if(t==0){r.adv=adv[0];r.diff=diff[0];r.bad=bad[0];r.total=total[0];r.cell=cell[0];results[blockIdx.x]=r;}
}
__global__ void rate_reduce_final(const Rate* input,int n,Rate* result) {
    __shared__ double totals[128],adv[128],diff[128],bad[128];__shared__ int cell[128];
    const int t=threadIdx.x;Rate local;
    for(int k=t;k<n;k+=blockDim.x)local=combine_rates(local,input[k]);
    totals[t]=local.total;adv[t]=local.adv;diff[t]=local.diff;bad[t]=local.bad;cell[t]=local.cell;__syncthreads();
    for(int stride=64;stride;stride/=2){if(t<stride){
        if(better_rate(totals[t+stride],cell[t+stride],totals[t],cell[t])){totals[t]=totals[t+stride];adv[t]=adv[t+stride];diff[t]=diff[t+stride];cell[t]=cell[t+stride];}
        bad[t]=fmax(bad[t],bad[t+stride]);}__syncthreads();}
    if(t==0){Rate r;r.total=totals[0];r.adv=adv[0];r.diff=diff[0];r.bad=bad[0];r.cell=cell[0];*result=r;}
}
#endif
