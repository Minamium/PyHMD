#pragma once
#include "pyhmd/physics/mhd.hpp"

HD inline double field1(const double* a,int i,int j,MeshView m,Parameters p) {
    return (a[m.vertex(i,j+1)]-a[m.vertex(i,j)])/m.area1(i,j,p);
}
HD inline double field2(const double* a,int i,int j,MeshView m,Parameters p) {
    double ar=m.area2(i,j,p);return fabs(ar)<1e-28?0.:-(a[m.vertex(i+1,j)]-a[m.vertex(i,j)])/ar;
}
