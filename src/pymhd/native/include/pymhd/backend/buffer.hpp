#pragma once
#include "pymhd/backend/execution.hpp"

template<class T> struct Buffer {
    T* p=nullptr; size_t n=0;
    Buffer()=default;
    explicit Buffer(size_t size) { allocate(size); }
    Buffer(const Buffer&)=delete; Buffer& operator=(const Buffer&)=delete;
    ~Buffer() {
#ifdef PYMHD_CUDA
        if(p) cudaFree(p);
#else
        delete[] p;
#endif
    }
    void allocate(size_t size) {
        if(p) throw std::runtime_error("buffer already allocated");
        n=size;
#ifdef PYMHD_CUDA
        cuda_check(cudaMalloc(reinterpret_cast<void**>(&p),sizeof(T)*n),"cudaMalloc");
#else
        p=new T[n];
#endif
    }
    void upload(const std::vector<T>& a) {
        if(a.size()!=n) throw std::runtime_error("buffer upload size mismatch");
#ifdef PYMHD_CUDA
        cuda_check(cudaMemcpy(p,a.data(),sizeof(T)*n,cudaMemcpyHostToDevice),"upload");
#else
        std::copy(a.begin(),a.end(),p);
#endif
    }
    std::vector<T> download() const {
        std::vector<T> a(n);
#ifdef PYMHD_CUDA
        cuda_check(cudaMemcpy(a.data(),p,sizeof(T)*n,cudaMemcpyDeviceToHost),"download");
#else
        std::copy(p,p+n,a.begin());
#endif
        return a;
    }
    void copy_from(const Buffer& a) {
        if(a.n!=n) throw std::runtime_error("buffer copy size mismatch");
#ifdef PYMHD_CUDA
        cuda_check(cudaMemcpyAsync(p,a.p,sizeof(T)*n,cudaMemcpyDeviceToDevice),"device copy");
#else
        std::copy(a.p,a.p+n,p);
#endif
    }
};
