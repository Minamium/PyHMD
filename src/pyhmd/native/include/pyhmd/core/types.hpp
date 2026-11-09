#pragma once
#include "pyhmd/backend/execution.hpp"

#define PYMHD_MODERN_THERMAL_SYNC 1
#define PYMHD_METRIC_CACHE 1
#define PYMHD_STAGGERED_CFL 1
constexpr int NV=10, NG=3;
constexpr double PI=3.14159265358979323846;
enum {RHO, M1, M2, M3, ENE, B1, B2, B3, ENTR, EINT};
enum Case {UNIFORM, ALFVEN, DIFFUSION, BRIO_WU, DISK};
enum Boundary {OUTFLOW, PERIODIC, REFLECT, AXIS, EQUATOR, PAPER_INNER, PAPER_OUTER};
struct State { double q[NV]{}; HD double& operator[](int i){return q[i];} HD double operator[](int i)const{return q[i];} };
struct Primitive { double q[NV]{}; HD double& operator[](int i){return q[i];} HD double operator[](int i)const{return q[i];} };
struct Current { double q[3]{}; HD double& operator[](int i){return q[i];} HD double operator[](int i)const{return q[i];} };
struct Weights { double w[4]{}; };
struct PpmGeometry;
struct Midplane { double va=0,cs=0,mu=0; };
struct EtaGeometry {int index=0;double fraction=0,height=0,profile=0;};
struct RadialGeometry {double volume=0,theta_area=0;};
struct PolarGeometry {double solid_angle=0,sine_face=0,sine_center=0,cosine_center=0;};
struct FaceSignal {double transverse_velocity=0,left=0,right=0,diffusion=0;};
// Diagnostic-only integrals. Face orientation is +coordinate; quantities are
// extensive per radian azimuth, including face area and accepted RK weights.
enum {TF_MASS,TF_MASS_ADV,TF_MASS_JUMP,TF_ENERGY,TF_ANGULAR,TF_ANG_ADV,TF_ANG_MAG,TF_ANG_JUMP,NTF};
struct TransportFace {double q[NTF]{};};
enum {TV_ADV,TV_JUMP,TV_RES,TV_ROUND,NTV};
struct TransportVertex {double q[NTV]{};};
struct HllTransport {double mass_adv=0,mass_jump=0,angular_adv=0,angular_mag=0,angular_jump=0;};
struct Parameters {
    int nx=600,ny=128,spherical=1,log_r=1,problem=DISK;
    int lower1=OUTFLOW,upper1=OUTFLOW,lower2=OUTFLOW,upper2=EQUATOR;
    int cooling=0,eta_model=0;
    // 0: total energy; 1: smooth entropy; 2: dual entropy; 3: dual internal energy.
    int entropy_switch=0,ct_energy_correction=0,max_step_retries=8,flat_reconstruction=0,reconstruction_cache=1;
    int metric_ppm=1;
    int initial_equilibrium=1; // 0: hydrostatic; 1: magnetohydrostatic disk.
    int retry_flat_same_dt=1;
    // Explicit reference-comparison option; strict rejection remains default.
    int negative_pressure_floor=0;
    // 0: inherited shock-neighbour veto; 1: Enzo-E modern local dual-energy
    // synchronization, meaningful only with entropy_switch==3.
    int internal_energy_sync=0;
    double dual_energy_eta=1e-3;
    int ct_emf_scheme=0; // 0: four-corner HLL; 1: face-Riemann HLL UCT.
    int cfl_staggered=1;
    int inner_flux_project=0; // Preserve physical normal CT flux; fill ghosts.
    int transport_diagnostics=0; // Passive accepted-step face/CT ledgers.
    double entropy_ratio=1e-3,entropy_shock_threshold=.05,entropy_contact_threshold=.1;
    double xlo=1,xhi=1500,ylo=0,yhi=PI/2,gamma=5.0/3,cfl=.4;
    double t_end=100,output_interval=10,gm=1,eta=.001;
    double rho_floor=1e-10,p_floor=1e-12,h=.1,corona=.001;
    double b0=.01414213562373095,mu0=.01,bend=.5,alpha_m=1.65,chi=.5,mu_ref=.01;
    double wave_amplitude=.1,glm_dummy=0;
    int max_steps=100000000,log_interval=100;
    const Midplane* midplane=nullptr;
    const Primitive *primitive_cache=nullptr,*edge_xl=nullptr,*edge_xr=nullptr,*edge_yl=nullptr,*edge_yr=nullptr;
    const int* shock_flags=nullptr;
    const Current* current_cache=nullptr;
    const EtaGeometry *eta_cell_geometry=nullptr,*eta_radial_geometry=nullptr,*eta_theta_geometry=nullptr,*eta_corner_geometry=nullptr;
};
struct MeshView {
    int nx,ny,tx,ty,start,global_nx;
    const double *xf,*yf,*xc,*yc,*qr,*qt;
    const Weights *wx,*wy;
    HD int cell(int i,int j)const{return j*tx+i;}
    HD int vertex(int i,int j)const{return j*(tx+1)+i;}
    HD double sine_center(int j)const{return polar_geometry?polar_geometry[j].sine_center:sin(yc[j]);}
    HD double cosine_center(int j)const{return polar_geometry?polar_geometry[j].cosine_center:cos(yc[j]);}
    HD double sine_face(int j)const{return polar_geometry?polar_geometry[j].sine_face:sin(yf[j]);}
    HD double radius(int i,int j,Parameters p)const{return p.spherical?xc[i]*sine_center(j):1.;}
    HD double volume(int i,int j,Parameters p)const {
        if(radial_geometry&&polar_geometry)return radial_geometry[i].volume*polar_geometry[j].solid_angle;
        return p.spherical?(pow(xf[i+1],3)-pow(xf[i],3))/3*(cos(yf[j])-cos(yf[j+1])):(xf[i+1]-xf[i])*(yf[j+1]-yf[j]);
    }
    HD double area1(int face,int j,Parameters p)const{return p.spherical?xf[face]*xf[face]*(polar_geometry?polar_geometry[j].solid_angle:cos(yf[j])-cos(yf[j+1])):yf[j+1]-yf[j];}
    HD double area2(int i,int face,Parameters p)const{return p.spherical?(radial_geometry?radial_geometry[i].theta_area:.5*(xf[i+1]*xf[i+1]-xf[i]*xf[i]))*sine_face(face):xf[i+1]-xf[i];}
    HD double edge_length(int i,int j,Parameters p)const{return p.spherical?xf[i]*sine_face(j):1.;}
    const Weights* wx_staggered=nullptr;
    const PpmGeometry *gx=nullptr,*gy=nullptr,*gx_staggered=nullptr;
    const RadialGeometry* radial_geometry=nullptr;
    const PolarGeometry* polar_geometry=nullptr;
};
