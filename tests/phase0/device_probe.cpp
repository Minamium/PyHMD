#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
int main() {
    int count=0;
    auto error=cudaGetDeviceCount(&count);
    if(error!=cudaSuccess){std::fprintf(stderr,"cudaGetDeviceCount: %s\n",cudaGetErrorString(error));return 2;}
    std::printf("{\"count\":%d,\"devices\":[",count);
    for(int device=0;device<count;++device){
        cudaDeviceProp p{};error=cudaGetDeviceProperties(&p,device);
        if(error!=cudaSuccess)return 2;
        if(device)std::printf(",");
        std::printf("{\"ordinal\":%d,\"name\":\"%s\",\"sm_count\":%d,\"major\":%d,\"minor\":%d,\"uuid_hex\":\"",device,p.name,p.multiProcessorCount,p.major,p.minor);
        for(unsigned char b:p.uuid.bytes)std::printf("%02x",b);
        std::printf("\"}");
    }
    std::printf("]}\n");return 0;
}
