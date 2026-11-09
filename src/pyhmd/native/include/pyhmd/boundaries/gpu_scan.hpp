#pragma once
#include "pyhmd/problems/builtin.hpp"

// One angular recurrence is an affine map y_j=a_j*y_(j-1)+b_j.
// Composition is associative; a prefix scan evaluates all angular positions
// concurrently. Radial ghost columns remain dependent; the physical normal
// flux is preserved unless the explicit legacy projection is selected.
struct PaperBoundaryAffine {double a,b;};
HD inline PaperBoundaryAffine paper_boundary_affine(const double* potential_a,int i,int j,MeshView m) {
    double t=.5*(m.yf[j]+m.yf[j-1]);
    double angle=(70*PI/180)/(1+exp(-(t-PI/4)/(PI/12)));
    double kval=m.xf[i]/((m.xf[i+1]-m.xf[i])*fmax(tan(angle),1e-10));
    double z=kval*(m.yf[j]-m.yf[j-1]),decay=exp(-z);
    double avg=z<1e-5?1-z/2+z*z/6-z*z*z/24:-expm1(-z)/z;
    return {decay,(avg-decay)*potential_a[m.vertex(i+1,j-1)]+(1-avg)*potential_a[m.vertex(i+1,j)]};
}
HD inline PaperBoundaryAffine compose_paper_boundary(PaperBoundaryAffine after,PaperBoundaryAffine before) {
    return {after.a*before.a,after.b+after.a*before.b};
}

#ifdef PYMHD_CUDA
__global__ void paper_inner_potential_scan(double* potential_a,MeshView m,int first_column) {
    constexpr int block=128;
    __shared__ double aa[block],bb[block],carry;
    int t=threadIdx.x;
    for(int i=first_column;i>=0;--i) {
        if(t==0){potential_a[m.vertex(i,NG)]=0;carry=0;}
        __syncthreads();
        for(int start=0;start<m.ny;start+=block) {
            int j=NG+start+t+1;
            PaperBoundaryAffine value=j<=NG+m.ny?paper_boundary_affine(potential_a,i,j,m):PaperBoundaryAffine{1,0};
            aa[t]=value.a;bb[t]=value.b;__syncthreads();
            for(int offset=1;offset<block;offset*=2) {
                value={aa[t],bb[t]};
                if(t>=offset)value=compose_paper_boundary(value,{aa[t-offset],bb[t-offset]});
                __syncthreads();aa[t]=value.a;bb[t]=value.b;__syncthreads();
            }
            double result=bb[t]+aa[t]*carry;
            if(j<=NG+m.ny)potential_a[m.vertex(i,j)]=result;
            __syncthreads();
            const int last=m.ny-start<block?m.ny-start-1:block-1;
            if(t==last)carry=result;
            __syncthreads();
        }
    }
}
inline void launch_paper_inner_potential(double* a,MeshView m,int first_column=NG-1) {
    paper_inner_potential_scan<<<1,128>>>(a,m,first_column);
    cuda_check(cudaGetLastError(),"paper inner potential scan");
}
#endif
