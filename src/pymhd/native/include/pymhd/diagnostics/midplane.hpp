#pragma once
#include "pymhd/numerics/thermal.hpp"
#include "pymhd/distributed/world.hpp"
#include <fstream>
#include <iomanip>

// Passive accepted-output sampling, independent of the face-flux ledger scheme.
// Psi is linearly interpolated at band endpoints; pressure is raw total-energy
// pressure, matching the existing transport-midplane diagnostic convention.
inline void write_midplane_sample(std::ofstream& stream,const World& world,
                                 const std::vector<State>& u,const std::vector<double>& psi,
                                 MeshView m,Parameters p,int step,double time) {
    double local[5]{},global[5]{};int j=NG+m.ny-1;
    for(int i=NG;i<NG+m.nx;++i) {
        double weight=std::max(0.,std::min(m.xf[i+1],1.5)-std::max(m.xf[i],1.1));
        double b=field2(psi.data(),i,j+1,m,p),pressure=raw_pressure(u[m.cell(i,j)],m.radius(i,j,p),p);
        local[0]+=weight*b*b/(2*pressure);local[1]+=weight*b*b;local[2]+=weight*pressure;local[3]+=weight;
        local[4]+=weight/(m.xf[i+1]-m.xf[i])*(psi[m.vertex(i+1,j+1)]-psi[m.vertex(i,j+1)]);
    }
#ifdef PYMHD_MPI
    MPI_Reduce(local,global,5,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
    std::copy(local,local+5,global);
#endif
    if(world.rank==0) {
        if(!(global[3]>0))throw std::runtime_error("empty midplane diagnostic band");
        stream<<std::setprecision(17)<<step<<','<<time;
        for(int k=0;k<3;++k)stream<<','<<global[k]/global[3];
        stream<<','<<global[3]<<','<<global[4]<<'\n';stream.flush();
        if(!stream)throw std::runtime_error("midplane metric output failed");
    }
}
