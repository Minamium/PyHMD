#pragma once
#include "pyhmd/backend/buffer.hpp"

struct World {
    int rank=0,size=1,local_rank=0;
    World(int& argc,char**& argv) {
#ifdef PYMHD_MPI
        MPI_Init(&argc,&argv); MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm local; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
        MPI_Comm_rank(local,&local_rank); MPI_Comm_free(&local);
#else
        (void)argc; (void)argv;
#endif
#ifdef PYMHD_CUDA
        int count=0; cuda_check(cudaGetDeviceCount(&count),"device count");
        if(count==0 || (count!=1 && local_rank>=count)) throw std::runtime_error("insufficient visible GPUs for local ranks");
        int device=count==1?0:local_rank; cuda_check(cudaSetDevice(device),"set device");
        cudaDeviceProp prop{}; cuda_check(cudaGetDeviceProperties(&prop,device),"device properties");
        if(const char* requested=std::getenv("PYMHD_GPU_BLOCK_SIZE")) {
            char* end=nullptr;long block=std::strtol(requested,&end,10);
            if(end==requested || *end!='\0' || block<32 || block>512 || block%32!=0 || block>prop.maxThreadsPerBlock)
                throw std::runtime_error("PYMHD_GPU_BLOCK_SIZE must be a warp multiple between 32 and 512 supported by the GPU");
            gpu_block_size()=static_cast<int>(block);
        }
#ifdef PYMHD_MPI
        // UUID uniqueness on each node; handles per-rank CUDA_VISIBLE_DEVICES masks.
        MPI_Comm node; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        int count_local=0; MPI_Comm_size(node,&count_local); std::vector<char> ids(16*count_local);
        MPI_Allgather(prop.uuid.bytes,16,MPI_CHAR,ids.data(),16,MPI_CHAR,node);
        for(int a=0;a<count_local;++a) for(int b=0;b<a;++b)
            if(std::memcmp(ids.data()+16*a,ids.data()+16*b,16)==0)
                throw std::runtime_error("two MPI ranks selected the same GPU UUID");
        MPI_Comm_free(&node);
#endif
        std::cerr<<"rank="<<rank<<" local_rank="<<local_rank<<" device="<<device<<" GPU="<<prop.name
                 <<" SMs="<<prop.multiProcessorCount<<" block="<<gpu_block_size()<<"\n";
#endif
    }
    ~World() {
#ifdef PYMHD_MPI
        MPI_Finalize();
#endif
    }
    double maximum(double v) const {
#ifdef PYMHD_MPI
        double r; MPI_Allreduce(&v,&r,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); return r;
#else
        return v;
#endif
    }
    void abort(int code) const {
#ifdef PYMHD_MPI
        MPI_Abort(MPI_COMM_WORLD,code);
#else
        (void)code;
#endif
    }
};
