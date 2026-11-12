#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef PYMHD_MPI
#include <mpi.h>
#endif
#ifdef PYMHD_CUDA
#include <cuda_runtime.h>
#include <typeinfo>
#define HD __host__ __device__
inline void cuda_check(cudaError_t e, const char* where) {
    if (e != cudaSuccess) throw std::runtime_error(std::string(where)+": "+cudaGetErrorString(e));
}
inline void synchronize() { cuda_check(cudaDeviceSynchronize(), "CUDA synchronize"); }
// Keep the launch geometry tunable: occupancy depends on the registers used by
// each kernel, and the largest block is not necessarily the fastest one.
inline int& gpu_block_size() { static int value=128;return value; }
template<class F> __global__ void parallel_kernel(int n, F f) {
    int k=blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n) f(k);
}
template<class F> void parallel_for(int n,F f) {
    if(n<=0)return;
    const int block=gpu_block_size();
    static const bool statistics=[](){const char* v=std::getenv("PYMHD_KERNEL_STATS");return v&&std::string(v)=="1";}();
    static bool reported=false;
    if(statistics&&!reported) {
        cudaFuncAttributes attributes{};int blocks=0,device=0;cudaDeviceProp properties{};
        cuda_check(cudaGetDevice(&device),"kernel stats device");
        cuda_check(cudaGetDeviceProperties(&properties,device),"kernel stats properties");
        cuda_check(cudaFuncGetAttributes(&attributes,parallel_kernel<F>),"kernel stats attributes");
        cuda_check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&blocks,parallel_kernel<F>,block,0),"kernel stats occupancy");
        std::cerr<<"kernel_stats name="<<typeid(F).name()<<" registers="<<attributes.numRegs
            <<" local_bytes="<<attributes.localSizeBytes<<" block="<<block<<" active_blocks_per_sm="<<blocks
            <<" theoretical_warp_occupancy="<<double(blocks*block)/properties.maxThreadsPerMultiProcessor<<'\n';
        reported=true;
    }
    parallel_kernel<<<(n+block-1)/block,block>>>(n,f);
    cuda_check(cudaGetLastError(),"CUDA kernel launch");
}
#else
#define HD
inline void synchronize() {}
template<class F> void parallel_for(int n,F f) { for(int k=0;k<n;++k) f(k); }
#endif
