// Independent moment, convergence and bounded-parabola checks. These test
// reconstruction accuracy, not the order/positivity of the complete MHD method.
#include "pymhd/setup.hpp"
#include "pymhd/metric_reconstruction.hpp"
#include <iomanip>
#include <random>

static void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
static long double jac(long double x,MetricMeasure m) {
    if(m==MetricMeasure::RadialVolume)return x*x;
    if(m==MetricMeasure::RadialArea)return std::fabs(x);
    if(m==MetricMeasure::PolarVolume)return std::fabs(std::sin(x));
    return 1;
}
// Independent composite Simpson integration, not the helper's Gauss rule.
template<class F> static long double average(double lo,double hi,MetricMeasure measure,F function) {
    constexpr int n=2048;long double numerator=0,denominator=0;
    for(int k=0;k<=n;++k) {
        const long double x=static_cast<long double>(lo)+(static_cast<long double>(hi)-lo)*k/n;
        const long double w=(k==0||k==n)?1:(k%2?4:2);
        const long double j=jac(x,measure);
        numerator+=w*j*function(x);denominator+=w*j;
    }
    return numerator/denominator;
}
static double interpolate(const std::vector<double>& faces,int face,MetricMeasure measure,
                          const std::vector<double>& values,bool metric) {
    std::vector<double> volume(faces.size());
    for(size_t i=0;i<faces.size();++i) {
        if(measure==MetricMeasure::RadialVolume)volume[i]=faces[i]*faces[i]*faces[i]/3;
        else if(measure==MetricMeasure::RadialArea)volume[i]=faces[i]*std::fabs(faces[i])/2;
        else if(measure==MetricMeasure::PolarVolume)volume[i]=setup_theta_volume(faces[i]);
        else volume[i]=faces[i];
    }
    Weights weights=metric?metric_face_weights(faces,face,measure):setup_face_weights(volume,face);
    double result=0;for(int k=0;k<4;++k)result+=weights.w[k]*values[face-2+k];return result;
}

int main() {
    try {
        std::cout<<std::setprecision(12);
        double max_polynomial_error=0,max_mean_error=0;
        // Uniform/irregular spherical radial meshes, reflected polar halos,
        // and small cells at large physical radius test different conditioning.
        for(MetricMeasure measure:{MetricMeasure::Cartesian,MetricMeasure::RadialArea,MetricMeasure::RadialVolume,MetricMeasure::PolarVolume})
        for(int family=0;family<3;++family) {
            std::vector<double> faces(17);
            for(int k=0;k<17;++k) {
                if(measure==MetricMeasure::PolarVolume)faces[k]=(k-3)*(.04/(family+1));
                else if(family==0)faces[k]=.03*k;
                else if(family==1)faces[k]=std::exp(.08*k);
                else faces[k]=1400+.001*k+.00001*k*k;
            }
            for(int face=2;face<=14;++face)for(int degree=0;degree<=3;++degree) {
                const long double origin=faces[face],scale=faces[face+2]-faces[face-2];
                std::vector<double> q(16);
                for(int c=0;c<16;++c)q[c]=static_cast<double>(average(faces[c],faces[c+1],measure,[&](long double x){return std::pow((x-origin)/scale,degree);}));
                double expected=degree==0?1.:0.;double actual=interpolate(faces,face,measure,q,true);
                max_polynomial_error=std::max(max_polynomial_error,std::fabs(actual-expected));
                // The independent absolute-coordinate quadrature loses a
                // few ulps/dx at r=1400 and dx~.001 on binary64 long double.
                const double tolerance=family==2&&measure!=MetricMeasure::PolarVolume?2e-10:2e-11;
                if(std::fabs(actual-expected)>=tolerance)std::cerr<<"measure="<<static_cast<int>(measure)<<" family="<<family<<" face="<<face<<" degree="<<degree<<" actual="<<std::setprecision(17)<<actual<<" expected="<<expected<<'\n';
                require(std::fabs(actual-expected)<tolerance,"physical-coordinate cubic moment reproduction failed");
            }
        }
        // Geometry-limit coefficients have known axis limits independent of
        // the general quadrature construction.
        auto radial=metric_ppm_geometry(0,1,MetricMeasure::RadialVolume);
        auto area=metric_ppm_geometry(0,1,MetricMeasure::RadialArea);
        auto polar=metric_ppm_geometry(0,1e-4,MetricMeasure::PolarVolume);
        require(std::fabs(radial.hplus-5)<1e-12&&std::fabs(radial.hminus-5./3)<1e-12,"incorrect spherical-origin parabola measure");
        require(std::fabs(area.hplus-4)<1e-12&&std::fabs(area.hminus-2)<1e-12,"incorrect face-area parabola measure");
        require(std::fabs(polar.hplus-4)<1e-8&&std::fabs(polar.hminus-2)<1e-8,"incorrect polar-axis parabola limit");

        std::mt19937_64 generator(14040537);std::uniform_real_distribution<double> random(-2,2);
        for(MetricMeasure measure:{MetricMeasure::Cartesian,MetricMeasure::RadialArea,MetricMeasure::RadialVolume,MetricMeasure::PolarVolume})
        for(auto cell:std::vector<std::pair<double,double>>{{0,.01},{.01,.023},{.3,.6},{1.4,1.5}}) {
            double lo=cell.first,hi=cell.second;auto geometry=metric_ppm_geometry(lo,hi,measure);
            const long double m1=average(lo,hi,measure,[&](long double x){return (x-lo)/(hi-lo);});
            const long double m12=average(lo,hi,measure,[&](long double x){long double y=(x-lo)/(hi-lo);return y*(1-y);});
            for(int sample=0;sample<80;++sample) {
                double qc=random(generator),left=random(generator),right=random(generator);
                const double old_left=left,old_right=right;
                metric_ppm_limit(qc,left,right,geometry);
                const double min=std::min({qc,old_left,old_right}),max=std::max({qc,old_left,old_right});
                const long double curve=-(geometry.hplus*(right-qc)+geometry.hminus*(left-qc));
                const long double mean=left+(right-left)*m1+curve*m12;
                max_mean_error=std::max(max_mean_error,std::fabs(static_cast<double>(mean)-qc));
                require(std::fabs(mean-qc)<2e-12,"limited parabola does not preserve weighted mean");
                long double previous=left;
                for(int n=0;n<=100;++n) {
                    const long double y=n/100.L,value=left+y*((right-left)+curve*(1-y));
                    require(value>=min-2e-12&&value<=max+2e-12,"metric PPM introduced an overshoot");
                    if(right>=left)require(value>=previous-2e-12,"limited parabola not monotone increasing");
                    else require(value<=previous+2e-12,"limited parabola not monotone decreasing");
                    previous=value;
                }
            }
        }
        // Smooth first non-axis interface: a geometric manufactured field.
        // The old reconstruction is polynomial in signed volume coordinate,
        // the improved one is polynomial in theta with |sin(theta)| moments.
        double previous_metric=0,previous_legacy=0;
        for(int n:{16,32,64,128}) {
            const double h=(PI/2)/n;std::vector<double> faces(9),values(8);
            for(int k=0;k<9;++k)faces[k]=(k-3)*h;
            for(int k=0;k<8;++k)values[k]=static_cast<double>(average(faces[k],faces[k+1],MetricMeasure::PolarVolume,[](long double x){return 2+.2L*std::sin(x);}));
            const int face=4;const double truth=2+.2*std::sin(h);
            // Reconstruct BOTH edges of the physical cell [0,h], bound by
            // adjacent averages, and apply the actual metric limiter.
            const int cell=face-1;
            double left=interpolate(faces,cell,MetricMeasure::PolarVolume,values,true);
            double right=interpolate(faces,cell+1,MetricMeasure::PolarVolume,values,true);
            left=std::max(std::min(values[cell-1],values[cell]),std::min(left,std::max(values[cell-1],values[cell])));
            right=std::max(std::min(values[cell],values[cell+1]),std::min(right,std::max(values[cell],values[cell+1])));
            metric_ppm_limit(values[cell],left,right,metric_ppm_geometry(faces[cell],faces[cell+1],MetricMeasure::PolarVolume));
            double metric=std::fabs(right-truth);
            double legacy=std::fabs(interpolate(faces,face,MetricMeasure::PolarVolume,values,false)-truth);
            std::cout<<"polar N="<<n<<" metric_error="<<metric<<" legacy_volume_error="<<legacy;
            if(previous_metric>1e-13)std::cout<<" metric_order="<<std::log(previous_metric/metric)/std::log(2.);
            if(previous_legacy>0)std::cout<<" legacy_order="<<std::log(previous_legacy/legacy)/std::log(2.);
            std::cout<<'\n';
            if(previous_metric>1e-12)require(metric<previous_metric/8,"metric reconstruction lost manufactured order near axis");
            require(metric<legacy*1e-2,"metric reconstruction did not improve polar first-interface accuracy");
            previous_metric=metric;previous_legacy=legacy;
        }
        // Distinguish face-area and volume measures: they are not interchangeable.
        std::vector<double> faces={0,1,2,3,4,5};std::vector<double> face_average(5);
        for(int c=0;c<5;++c)face_average[c]=static_cast<double>(average(faces[c],faces[c+1],MetricMeasure::RadialArea,[](long double r){return r;}));
        require(std::fabs(interpolate(faces,2,MetricMeasure::RadialArea,face_average,true)-2)<1e-12,"staggered radial-area linear field not reproduced");
        require(std::fabs(interpolate(faces,2,MetricMeasure::RadialVolume,face_average,true)-2)>1e-3,"test failed to distinguish incorrect staggered measure");
        std::cout<<"PASS metric reconstruction: maximum cubic error="<<max_polynomial_error<<", maximum weighted-mean error="<<max_mean_error<<'\n';
    }catch(const std::exception& error){std::cerr<<"FAIL metric reconstruction: "<<error.what()<<'\n';return 1;}
}
