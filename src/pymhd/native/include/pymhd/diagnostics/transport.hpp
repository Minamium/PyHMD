#pragma once
#include "pymhd/runtime/solver.hpp"
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>

// Host-side output of passive accepted-step ledgers. No live boundary fill or
// operator evaluation is performed here. Region cuts use a 600x128 reference
// grid, making their physical faces identical on the nested 1200x256 grid.
class TransportOutput {
    struct Region {std::string name;int i0,i1,j0,j1;};
    enum {M,E,L,M0,E0,L0,RADIAL,POLAR=RADIAL+NTF,FM=POLAR+NTF,FE,FP,CT,INTERNAL,ENTROPY,BOUNDARY,COOLING,VOLUME,NS};
    const World& world;
    MeshView m;
    Parameters p;
    std::filesystem::path directory;
    std::vector<Region> regions;
    std::ofstream stocks,boundary,midplane;
    static constexpr const char* names[NTF]={"mass","mass_adv","mass_hll_jump","energy","angular","angular_adv","angular_maxwell","angular_hll_jump"};
    int radial_cut(double radius)const {
        int reference=static_cast<int>(std::llround(std::log(radius/p.xlo)/std::log(p.xhi/p.xlo)*600));
        return std::max(0,std::min(p.nx,static_cast<int>(std::llround(reference*p.nx/600.))));
    }
    int theta_cut()const {
        int reference=static_cast<int>(std::llround((PI/2-2*p.h-p.ylo)/(p.yhi-p.ylo)*128));
        return std::max(0,std::min(p.ny,static_cast<int>(std::llround(reference*p.ny/128.))));
    }
    double radial_face(int i)const{return p.xlo*std::exp(std::log(p.xhi/p.xlo)*i/p.nx);}
    double theta_face(int j)const{return p.ylo+(p.yhi-p.ylo)*j/p.ny;}
public:
    TransportOutput(const std::string& output,const World& w,MeshView host,Parameters params):world(w),m(host),p(params),directory(output) {
        const int cut=theta_cut();
        regions={{"mu_disk",radial_cut(1.1),radial_cut(1.5),cut,p.ny},
                 {"disk_r10",0,radial_cut(10),cut,p.ny},
                 {"disk_r100",0,radial_cut(100),cut,p.ny},
                 {"corona_r10",0,radial_cut(10),0,cut},
                 {"domain",0,p.nx,0,p.ny}};
        if(world.rank)return;
        std::filesystem::create_directories(directory);
        stocks.open(directory/"transport_regions.csv");
        stocks<<"step,time,region,r_inner,r_outer,theta_lower,theta_upper,mass,energy,angular,mass_initial,energy_initial,angular_initial";
        for(const char* direction:{"radial","polar"})for(const char* name:names)stocks<<','<<direction<<'_'<<name<<"_out";
        stocks<<",floor_mass,floor_energy,floor_potential,ct_energy,internal_energy,entropy_energy,boundary_energy,cooling,volume,mass_residual,energy_residual,angular_residual\n";
        boundary.open(directory/"transport_inner_boundary.csv");
        boundary<<"step,time,j,theta,is_source_disk,active_rho,active_vr,ghost_rho,ghost_vr";
        for(const char* name:names)boundary<<",cumulative_"<<name<<"_out";
        boundary<<'\n';
        midplane.open(directory/"transport_midplane.csv");
        midplane<<"step,time,mu_r1p1_r1p5,Btheta_squared_mean,pressure_mean,density_mean,va_mean,eta_p_mean,radial_weight,psi_change_closure_max,emf_operator_round_sum_abs\n";
        std::ofstream meta(directory/"transport_method.json");
        meta<<"{\"schema\":1,\"normalization\":\"per radian azimuth, one simulated hemisphere\","
            <<"\"origin_time\":0,\"origin_step\":0,\"rk_weights\":[0.16666666666666667,0.16666666666666667,0.66666666666666667],"
            <<"\"mass_adv\":\"HLL wave-weighted left/right advective flux\","
            <<"\"ct_jump\":\"explicit UCT jump contribution; not the full discretization error\","
            <<"\"psi_change_closure\":\"actual Psi change minus sum of integrated operator terms; not absorbed into operator roundoff\","
            <<"\"region_boundary\":\"mesh-aligned cuts anchored to nearest faces of reference 600x128 grid\","
            <<"\"shared_vertices\":\"per-rank files retain duplicates; lowest radial rank owns the shared vertex\","
            <<"\"snapshot_faces\":\"per-rank integrated numerical flux files retain duplicate MPI interfaces\","
            <<"\"source_disk_boundary\":\"theta cell centre > pi/2 - 2 epsilon; distinct from mesh-aligned control cut\","
            <<"\"rates\":\"time differences of cumulative accepted fluxes divided by elapsed physical time\"}\n";
    }
    void write(Solver& s,const std::vector<State>& u,const std::vector<double>& psi,
               const std::vector<State>& ghosts,int step,double time,bool snapshot) {
        s.ensure_transport_origin(time,step);
        auto fr=s.transport_radial.download(),ft=s.transport_polar.download();
        auto pv=s.transport_psi.download();auto u0=s.transport_origin_state.download();auto psi0=s.transport_origin_psi.download();
        auto fm=s.floor_mass.download(),fe=s.floor_energy.download(),be=s.boundary_energy.download();
        auto thermal=s.thermal.download();auto floor_potential=s.flux_accounting.download();
        std::vector<double> values(regions.size()*NS,0),combined(values.size(),0);
        for(std::size_t n=0;n<regions.size();++n) {
            const auto& region=regions[n];double* q=values.data()+n*NS;
            for(int j=region.j0;j<region.j1;++j)for(int i=0;i<m.nx;++i) {
                int global=m.start+i;if(global<region.i0||global>=region.i1)continue;
                int id=m.cell(i+NG,j+NG),k=j*m.nx+i;double vol=m.volume(i+NG,j+NG,p),phi=potential(m.xc[i+NG],p);
                q[M]+=u[id][RHO]*vol;q[E]+=(u[id][ENE]+phi*u[id][RHO])*vol;q[L]+=u[id][M3]*vol;
                q[M0]+=u0[id][RHO]*vol;q[E0]+=(u0[id][ENE]+phi*u0[id][RHO])*vol;q[L0]+=u0[id][M3]*vol;
                for(int c=0;c<NTF;++c) {
                    if(global==region.i0)q[RADIAL+c]-=fr[j*(m.nx+1)+i].q[c];
                    if(global==region.i1-1)q[RADIAL+c]+=fr[j*(m.nx+1)+i+1].q[c];
                    if(j==region.j0)q[POLAR+c]-=ft[j*m.nx+i].q[c];
                    if(j==region.j1-1)q[POLAR+c]+=ft[(j+1)*m.nx+i].q[c];
                }
                q[FM]+=fm[k];q[FE]+=fe[k];q[FP]+=floor_potential[k].floor_potential;
                q[CT]+=thermal[k].ct_energy;q[INTERNAL]+=thermal[k].internal_energy;q[ENTROPY]+=thermal[k].entropy_energy;
                q[BOUNDARY]+=be[k];q[COOLING]+=thermal[k].compatible_cooling;q[VOLUME]+=vol;
            }
        }
        double profile[7]{},global_profile[7]{},closure=0,global_closure=0,operator_round=0,global_round=0;
        for(int i=NG;i<NG+m.nx;++i) {
            const int j=NG+m.ny-1;double weight=std::max(0.,std::min(m.xf[i+1],1.5)-std::max(m.xf[i],1.1));
            double btheta=field2(psi.data(),i,j+1,m,p),pressure=raw_pressure(u[m.cell(i,j)],m.radius(i,j,p),p),rho=u[m.cell(i,j)][RHO];
            profile[0]+=weight*btheta*btheta/(2*pressure);profile[1]+=weight*btheta*btheta;profile[2]+=weight*pressure;
            profile[3]+=weight*rho;profile[4]+=weight*std::abs(btheta)/std::sqrt(rho);
            double va=std::abs(btheta)/std::sqrt(rho),mu=btheta*btheta/(2*pressure),height=2*p.h*m.xc[i];
            double eta=p.eta_model==0?p.eta:p.alpha_m*va*height;
            if(p.eta_model==2)eta=p.alpha_m*std::sqrt(2*p.mu_ref)*sq(mu/p.mu_ref)*std::sqrt(pressure/rho)*height;
            profile[5]+=weight*eta;profile[6]+=weight;
        }
        for(int j=NG;j<=NG+m.ny;++j)for(int i=NG;i<=NG+m.nx;++i) {
            if(world.rank>0&&i==NG)continue;
            int v=m.vertex(i,j);double sum=0;for(double x:pv[v].q)sum+=x;
            closure=std::max(closure,std::abs((psi[v]-psi0[v])-sum));operator_round+=std::abs(pv[v].q[TV_ROUND]);
        }
#ifdef PYMHD_MPI
        MPI_Reduce(values.data(),combined.data(),static_cast<int>(values.size()),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(profile,global_profile,7,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&closure,&global_closure,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        MPI_Reduce(&operator_round,&global_round,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
        combined=values;std::copy(profile,profile+7,global_profile);global_closure=closure;global_round=operator_round;
#endif
        if(world.rank==0) {
            stocks<<std::setprecision(17);boundary<<std::setprecision(17);midplane<<std::setprecision(17);
            for(std::size_t n=0;n<regions.size();++n) {
                const auto& r=regions[n];const double* q=combined.data()+n*NS;
                stocks<<step<<','<<time<<','<<r.name<<','<<radial_face(r.i0)<<','<<radial_face(r.i1)<<','<<theta_face(r.j0)<<','<<theta_face(r.j1);
                for(int c=0;c<NS;++c)stocks<<','<<q[c];
                stocks<<','<<(q[M]-q[M0]+q[RADIAL+TF_MASS]+q[POLAR+TF_MASS]-q[FM])
                      <<','<<(q[E]-q[E0]+q[RADIAL+TF_ENERGY]+q[POLAR+TF_ENERGY]+q[COOLING]-q[FE]-q[FP]-q[CT]-q[INTERNAL]-q[ENTROPY]-q[BOUNDARY])
                      <<','<<(q[L]-q[L0]+q[RADIAL+TF_ANGULAR]+q[POLAR+TF_ANGULAR])<<'\n';
            }
            for(int j=0;j<m.ny;++j) {
                auto active=get_primitive(u.data(),NG,j+NG,m,p);
                auto ghost=get_primitive(ghosts.data(),NG-1,j+NG,m,p);
                boundary<<step<<','<<time<<','<<j<<','<<m.yc[j+NG]<<','<<(m.yc[j+NG]>PI/2-2*p.h?1:0)
                        <<','<<active[RHO]<<','<<active[M1]<<','<<ghost[RHO]<<','<<ghost[M1];
                for(double x:fr[j*(m.nx+1)].q)boundary<<','<<-x;
                boundary<<'\n';
            }
            midplane<<step<<','<<time;
            for(int c=0;c<6;++c)midplane<<','<<(global_profile[6]>0?global_profile[c]/global_profile[6]:0);
            midplane<<','<<global_profile[6]<<','<<global_closure<<','<<global_round<<'\n';
            stocks.flush();boundary.flush();midplane.flush();
            if(!stocks||!boundary||!midplane)throw std::runtime_error("transport CSV write failed");
        }
        if(snapshot) {
            std::ostringstream suffix;suffix<<"_step"<<std::setw(10)<<std::setfill('0')<<step<<"_rank"<<std::setw(6)<<world.rank<<".csv";
            const auto folder=directory/"transport_snapshots";std::filesystem::create_directories(folder);
            std::ofstream vertices(folder/("psi"+suffix.str()));vertices<<std::setprecision(17);
            vertices<<"# time="<<time<<" step="<<step<<" origin_time="<<s.transport_origin_time<<"\ni,j,r_face,theta_face,owner,psi,psi_initial,cumulative_advection,cumulative_uct_jump,cumulative_resistive,cumulative_operator_round,psi_change_residual\n";
            for(int j=NG;j<=NG+m.ny;++j)for(int i=NG;i<=NG+m.nx;++i) {
                int v=m.vertex(i,j);double sum=0;for(double x:pv[v].q)sum+=x;
                vertices<<m.start+i-NG<<','<<j-NG<<','<<m.xf[i]<<','<<m.yf[j]<<','<<(world.rank==0||i>NG?1:0)<<','<<psi[v]<<','<<psi0[v];
                for(double x:pv[v].q)vertices<<','<<x;
                vertices<<','<<((psi[v]-psi0[v])-sum)<<'\n';
            }
            std::ofstream faces(folder/("faces"+suffix.str()));faces<<std::setprecision(17);
            faces<<"# time="<<time<<" step="<<step<<"\ndir,i,j,owner";
            for(const char* name:names)faces<<",cumulative_"<<name<<"_plus_coordinate";
            faces<<'\n';
            for(int dir=0;dir<2;++dir) {
                const auto& f=dir==0?fr:ft;int ni=dir==0?m.nx+1:m.nx,nj=dir==0?m.ny:m.ny+1;
                for(int j=0;j<nj;++j)for(int i=0;i<ni;++i) {
                    faces<<dir<<','<<m.start+i<<','<<j<<','<<(dir!=0||world.rank==0||i>0?1:0);
                    for(double x:f[j*ni+i].q)faces<<','<<x;
                    faces<<'\n';
                }
            }
            vertices.flush();faces.flush();if(!vertices||!faces)throw std::runtime_error("transport snapshot write failed");
        }
    }
};
