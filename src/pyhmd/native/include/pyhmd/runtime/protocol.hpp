#pragma once
#include "pyhmd/config/native.hpp"
#include <iostream>
#ifndef PYMHD_SOURCE_FINGERPRINT
#define PYMHD_SOURCE_FINGERPRINT "unversioned-manual-build"
#endif

namespace pyhmd {
inline void describe() {
    std::cout<<"{\"protocol_version\":1,\"application\":\"PyHMD\",\"dtype\":\"float64\",\"dimensions\":2,"
             <<"\"model\":\"resistive_mhd\",\"scheme\":\"ppm_hll_ssprk3\",\"geometries\":[\"cartesian\",\"spherical\"],"
             <<"\"source_fingerprint\":\""<<PYMHD_SOURCE_FINGERPRINT<<"\",\"backend\":\""
#ifdef PYMHD_CUDA
             <<"cuda"
#else
             <<"cpu"
#endif
             <<"\",\"mpi\":"
#ifdef PYMHD_MPI
             <<"true"
#else
             <<"false"
#endif
             <<",\"fp_mode\":\""
#ifdef PYHMD_STRICT_FP
             <<"strict"
#else
             <<"default"
#endif
             <<"\",\"problems\":[\"thin_disk\",\"uniform\",\"alfven\",\"diffusion\",\"brio_wu\",\"oblique_alfven\",\"linear_wave\",\"field_loop\",\"rj2a\",\"orszag_tang\",\"blast_low_beta\"],\"wave_modes\":[\"fast\",\"alfven\",\"slow\",\"entropy\"]}\n";
}
inline int validate_configuration(int argc,char** argv) {
    try {
        auto setup=parse_setup(argc,argv);
        std::cout<<"{\"valid\":true,\"nx1\":"<<setup.p.nx<<",\"nx2\":"<<setup.p.ny<<"}\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
}
