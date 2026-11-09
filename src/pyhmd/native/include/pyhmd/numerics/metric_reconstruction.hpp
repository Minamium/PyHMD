#pragma once
#include "pyhmd/physics/mhd.hpp"
#include <limits>

// Modern, explicitly selected curvilinear reconstruction, not a claim of
// PLUTO 4.0 source equivalence. See Mignone (2014), arXiv:1404.0537,
// Eqs. (16), (20), (21), (46), (47).
// A scalar spherical radial cell average has J=r^2. Btheta stored on a
// constant-theta face instead has radial face-area measure J=r.
enum class MetricMeasure { Cartesian, RadialArea, RadialVolume, PolarVolume };
struct PpmGeometry { double hplus=3.,hminus=3.; };

namespace metric_reconstruction_detail {
inline long double jacobian(long double x,MetricMeasure measure) {
    switch(measure) {
        case MetricMeasure::Cartesian:return 1.;
        case MetricMeasure::RadialArea:return std::fabs(x);
        case MetricMeasure::RadialVolume:return x*x;
        // Reflect positive cell measures across polar-axis ghost zones.
        case MetricMeasure::PolarVolume:return std::fabs(std::sin(x));
    }
    throw std::runtime_error("unknown reconstruction measure");
}
// Eight-point Gauss-Legendre quadrature. Radial polynomial moments of degree
// at most five are exact; trigonometric moments are resolved to roundoff at
// the cell widths used here. All powers use scaled, local coordinates.
template<class F> inline long double integral(double lower,double upper,MetricMeasure measure,F function) {
    if(!std::isfinite(lower)||!std::isfinite(upper)||!(upper>lower))
        throw std::runtime_error("invalid metric reconstruction cell");
    static constexpr long double node[4]={
        .183434642495649804939476142360184L,.525532409916328985817739049189246L,
        .796666477413626739591553936475831L,.960289856497536231683560868569473L};
    static constexpr long double weight[4]={
        .362683783378361982965150449277195L,.313706645877887287337962201986601L,
        .222381034453374470544355994426240L,.101228536290376259152531354309962L};
    const long double width=static_cast<long double>(upper)-lower;
    const long double half=width/2;
    long double sum=0;
    for(int n=0;n<4;++n)for(int sign:{-1,1}) {
        const long double y=(1+sign*node[n])/2;
        const long double x=lower+width*y;
        sum+=weight[n]*jacobian(x,measure)*function(x,y);
    }
    return half*sum;
}
inline long double power(long double x,int degree) {
    long double result=1;for(int n=0;n<degree;++n)result*=x;return result;
}
} // namespace metric_reconstruction_detail

// Four-cell cubic polynomial in PHYSICAL x, constrained by weighted averages.
// The return stencil is [face-2,face-1,face,face+1], as for legacy Weights.
inline Weights metric_face_weights(const std::vector<double>& faces,int face,MetricMeasure measure) {
    Weights result;
    const int count=static_cast<int>(faces.size())-1;
    if(face<2||face>count-2)return result;
    using namespace metric_reconstruction_detail;
    const long double origin=faces[face],scale=static_cast<long double>(faces[face+2])-faces[face-2];
    if(!(scale>0)||!std::isfinite(scale))throw std::runtime_error("invalid metric interpolation stencil");
    long double matrix[4][5]{};
    for(int cell=0;cell<4;++cell) {
        const int n=face-2+cell;
        const long double volume=integral(faces[n],faces[n+1],measure,[](long double,long double){return 1.L;});
        if(!(volume>0)||!std::isfinite(volume))throw std::runtime_error("nonpositive metric cell measure");
        matrix[0][cell]=1;
        const long double offset=(static_cast<long double>(faces[n])-origin)/scale;
        const long double width=(static_cast<long double>(faces[n+1])-faces[n])/scale;
        for(int degree=1;degree<4;++degree)
            matrix[degree][cell]=integral(faces[n],faces[n+1],measure,[&](long double,long double y){return power(offset+width*y,degree);})/volume;
    }
    matrix[0][4]=1;
    for(int pivot=0;pivot<4;++pivot) {
        int row=pivot;
        for(int other=pivot+1;other<4;++other)
            if(std::fabs(matrix[other][pivot])>std::fabs(matrix[row][pivot]))row=other;
        if(std::fabs(matrix[row][pivot])<1e-24L)throw std::runtime_error("singular metric moment stencil");
        if(row!=pivot)for(int col=pivot;col<5;++col)std::swap(matrix[row][col],matrix[pivot][col]);
        const long double divisor=matrix[pivot][pivot];
        for(int col=pivot;col<5;++col)matrix[pivot][col]/=divisor;
        for(int other=0;other<4;++other)if(other!=pivot) {
            const long double factor=matrix[other][pivot];
            for(int col=pivot;col<5;++col)matrix[other][col]-=factor*matrix[pivot][col];
        }
    }
    for(int cell=0;cell<4;++cell)result.w[cell]=static_cast<double>(matrix[cell][4]);
    return result;
}

inline PpmGeometry metric_ppm_geometry(double lower,double upper,MetricMeasure measure) {
    if(measure==MetricMeasure::Cartesian) {
        if(!std::isfinite(lower)||!std::isfinite(upper)||!(upper>lower))throw std::runtime_error("invalid PPM cell");
        return {};
    }
    using namespace metric_reconstruction_detail;
    const long double m1=integral(lower,upper,measure,[](long double,long double y){return y;});
    const long double m0=integral(lower,upper,measure,[](long double,long double){return 1.L;});
    const long double m12=integral(lower,upper,measure,[](long double,long double y){return y*(1-y);});
    if(!(m0>0)||!(m12>0))throw std::runtime_error("invalid PPM weighted moments");
    PpmGeometry g{static_cast<double>(m1/m12),static_cast<double>((m0-m1)/m12)};
    if(!std::isfinite(g.hplus)||!std::isfinite(g.hminus)||!(g.hplus>1)||!(g.hminus>1))
        throw std::runtime_error("invalid PPM curvature coefficients");
    return g;
}

// Call after bounding each interface value by its two adjacent averages.
// This preserves the given weighted mean in the implicit cell parabola and
// removes interior extrema. As with standard PPM it flattens smooth extrema;
// no claim of high order at extrema or positivity of the full MHD step follows.
HD inline void metric_ppm_limit(double average,double& left,double& right,PpmGeometry geometry) {
    const double dm=left-average,dp=right-average;
    if((dm>=0&&dp>=0)||(dm<=0&&dp<=0)) {left=right=average;return;}
    const double cap_plus=(geometry.hminus+1)/(geometry.hplus-1);
    const double cap_minus=(geometry.hplus+1)/(geometry.hminus-1);
    if(fabs(dp)>=cap_plus*fabs(dm))right=average-cap_plus*dm;
    else if(fabs(dm)>=cap_minus*fabs(dp))left=average-cap_minus*dp;
}
