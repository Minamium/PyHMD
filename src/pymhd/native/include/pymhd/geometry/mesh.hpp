#pragma once
#include "pymhd/config/native.hpp"
#include "pymhd/numerics/metric_reconstruction.hpp"
#include "pymhd/distributed/world.hpp"

inline double setup_theta_volume(double theta) {
    if(theta<0)return std::cos(theta)-2;
    if(theta>PI)return 2+std::cos(theta);
    return -std::cos(theta);
}

inline Weights setup_face_weights(const std::vector<double>& coordinate,int face) {
    Weights weights;
    const int count=static_cast<int>(coordinate.size())-1;
    if(face<2 || face>count-2)return weights;
    const double origin=coordinate[face];
    const double scale=coordinate[face+2]-coordinate[face-2];
    if(!(scale>0))throw std::runtime_error("nonmonotone volume-coordinate interpolation stencil");
    double matrix[4][5]{};
    for(int degree=0;degree<4;++degree) {
        for(int cell=0;cell<4;++cell) {
            const int index=face-2+cell;
            const double a=(coordinate[index]-origin)/scale,b=(coordinate[index+1]-origin)/scale;
            if(!(b>a))throw std::runtime_error("degenerate PPM cell-average stencil");
            matrix[degree][cell]=(std::pow(b,degree+1)-std::pow(a,degree+1))/((degree+1)*(b-a));
        }
        matrix[degree][4]=degree==0?1:0;
    }
    for(int pivot=0;pivot<4;++pivot) {
        int row=pivot;
        for(int candidate=pivot+1;candidate<4;++candidate)if(std::abs(matrix[candidate][pivot])>std::abs(matrix[row][pivot]))row=candidate;
        if(std::abs(matrix[row][pivot])<1e-15)throw std::runtime_error("singular PPM moment matrix");
        if(row!=pivot)for(int column=pivot;column<5;++column)std::swap(matrix[row][column],matrix[pivot][column]);
        const double value=matrix[pivot][pivot];
        for(int column=pivot;column<5;++column)matrix[pivot][column]/=value;
        for(int other=0;other<4;++other)if(other!=pivot) {
            const double factor=matrix[other][pivot];
            for(int column=pivot;column<5;++column)matrix[other][column]-=factor*matrix[pivot][column];
        }
    }
    for(int cell=0;cell<4;++cell)weights.w[cell]=matrix[cell][4];
    return weights;
}

struct Mesh {
    Buffer<double> xf,yf,xc,yc,qr,qt;
    Buffer<Weights> wx,wy,wx_staggered;
    Buffer<PpmGeometry> gx,gy,gx_staggered;
    Buffer<RadialGeometry> radial_geometry;
    Buffer<PolarGeometry> polar_geometry;
    MeshView view{};
    std::vector<double> global_xfaces,global_yfaces;
    Mesh(Parameters p,const World& world) {
        const int base=p.nx/world.size,remainder=p.nx%world.size;
        const int nx=base+(world.rank<remainder?1:0),start=world.rank*base+std::min(world.rank,remainder);
        if(nx<NG)throw std::runtime_error("each MPI rank must own at least 3 radial cells; reduce ranks or increase nx1");
        const int tx=nx+2*NG,ty=p.ny+2*NG;
        const auto radial_face=[&](int index) {
            return p.log_r?std::exp(std::log(p.xlo)+static_cast<double>(index)/p.nx*std::log(p.xhi/p.xlo)):
                           p.xlo+(p.xhi-p.xlo)*static_cast<double>(index)/p.nx;
        };
        const auto polar_face=[&](int index) {return p.ylo+(p.yhi-p.ylo)*static_cast<double>(index)/p.ny;};
        global_xfaces.resize(p.nx+1);global_yfaces.resize(p.ny+1);
        for(int i=0;i<=p.nx;++i)global_xfaces[i]=radial_face(i);
        for(int j=0;j<=p.ny;++j)global_yfaces[j]=polar_face(j);
        std::vector<double> host_xf(tx+1),host_yf(ty+1),host_xc(tx),host_yc(ty),host_qr(tx),host_qt(ty),volume_x(tx+1),volume_y(ty+1);
        for(int i=0;i<=tx;++i) {
            host_xf[i]=radial_face(start+i-NG);
            if(p.spherical && host_xf[i]<=0)throw std::runtime_error("radial halo crosses r=0; use log radial spacing or a larger inner radius");
            volume_x[i]=p.spherical?std::pow(host_xf[i],3)/3:host_xf[i];
        }
        for(int j=0;j<=ty;++j) {
            host_yf[j]=polar_face(j-NG);
            volume_y[j]=p.spherical?setup_theta_volume(host_yf[j]):host_yf[j];
        }
        for(int i=0;i<tx;++i) {
            host_xc[i]=p.spherical?std::cbrt(.5*(std::pow(host_xf[i],3)+std::pow(host_xf[i+1],3))):.5*(host_xf[i]+host_xf[i+1]);
            host_qr[i]=p.spherical?.5*(host_xf[i+1]*host_xf[i+1]-host_xf[i]*host_xf[i])/(volume_x[i+1]-volume_x[i]):0;
        }
        for(int j=0;j<ty;++j) {
            host_yc[j]=.5*(host_yf[j]+host_yf[j+1]);
            const double denominator=std::cos(host_yf[j])-std::cos(host_yf[j+1]);
            host_qt[j]=p.spherical && std::abs(denominator)>1e-30?(std::sin(host_yf[j+1])-std::sin(host_yf[j]))/denominator:0;
        }
        std::vector<Weights> host_wx(tx+1),host_wy(ty+1),host_ws(tx+1);
        std::vector<PpmGeometry> host_gx(tx),host_gy(ty),host_gs(tx);
        const auto radial=p.spherical?MetricMeasure::RadialVolume:MetricMeasure::Cartesian;
        const auto staggered=p.spherical?MetricMeasure::RadialArea:MetricMeasure::Cartesian;
        const auto polar=p.spherical?MetricMeasure::PolarVolume:MetricMeasure::Cartesian;
        for(int i=0;i<=tx;++i) {
            host_wx[i]=p.metric_ppm?metric_face_weights(host_xf,i,radial):setup_face_weights(volume_x,i);
            host_ws[i]=p.metric_ppm?metric_face_weights(host_xf,i,staggered):host_wx[i];
        }
        for(int j=0;j<=ty;++j)host_wy[j]=p.metric_ppm?metric_face_weights(host_yf,j,polar):setup_face_weights(volume_y,j);
        for(int i=0;i<tx;++i) {
            host_gx[i]=p.metric_ppm?metric_ppm_geometry(host_xf[i],host_xf[i+1],radial):PpmGeometry{};
            host_gs[i]=p.metric_ppm?metric_ppm_geometry(host_xf[i],host_xf[i+1],staggered):PpmGeometry{};
        }
        for(int j=0;j<ty;++j)host_gy[j]=p.metric_ppm?metric_ppm_geometry(host_yf[j],host_yf[j+1],polar):PpmGeometry{};
        xf.allocate(host_xf.size());yf.allocate(host_yf.size());xc.allocate(host_xc.size());yc.allocate(host_yc.size());
        qr.allocate(host_qr.size());qt.allocate(host_qt.size());wx.allocate(host_wx.size());wy.allocate(host_wy.size());
        wx_staggered.allocate(host_ws.size());gx.allocate(host_gx.size());gy.allocate(host_gy.size());gx_staggered.allocate(host_gs.size());
        xf.upload(host_xf);yf.upload(host_yf);xc.upload(host_xc);yc.upload(host_yc);
        qr.upload(host_qr);qt.upload(host_qt);wx.upload(host_wx);wy.upload(host_wy);
        wx_staggered.upload(host_ws);gx.upload(host_gx);gy.upload(host_gy);gx_staggered.upload(host_gs);
        view={nx,p.ny,tx,ty,start,p.nx,xf.p,yf.p,xc.p,yc.p,qr.p,qt.p,wx.p,wy.p,wx_staggered.p,gx.p,gy.p,gx_staggered.p};
        // Immutable separable metrics avoid repeated FP64 powers and trig in
        // face/corner/RK cells without changing the finite-volume quadrature.
        const char* use_metrics=std::getenv("PYMHD_METRIC_CACHE");
        if(!use_metrics||std::string(use_metrics)!="0") {
            std::vector<RadialGeometry> rg(tx);std::vector<PolarGeometry> pg(ty+1);
            for(int i=0;i<tx;++i){rg[i].volume=p.spherical?(std::pow(host_xf[i+1],3)-std::pow(host_xf[i],3))/3:host_xf[i+1]-host_xf[i];rg[i].theta_area=.5*(host_xf[i+1]*host_xf[i+1]-host_xf[i]*host_xf[i]);}
            for(int j=0;j<=ty;++j){pg[j].sine_face=std::sin(host_yf[j]);if(j<ty){pg[j].solid_angle=p.spherical?std::cos(host_yf[j])-std::cos(host_yf[j+1]):host_yf[j+1]-host_yf[j];pg[j].sine_center=std::sin(host_yc[j]);pg[j].cosine_center=std::cos(host_yc[j]);}}
            radial_geometry.allocate(rg.size());polar_geometry.allocate(pg.size());radial_geometry.upload(rg);polar_geometry.upload(pg);
            view.radial_geometry=radial_geometry.p;view.polar_geometry=polar_geometry.p;
        }
    }
};
