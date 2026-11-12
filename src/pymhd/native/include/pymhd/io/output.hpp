#pragma once
#include "pymhd/config/native.hpp"
#include "pymhd/problems/initialize.hpp"
#include "pymhd/io/checkpoint.hpp"
#include "pymhd/diagnostics/transport.hpp"
#include "pymhd/diagnostics/midplane.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>

namespace fs=std::filesystem;
static std::string case_name(int c) {switch(c){case DISK:return "thin_disk";case UNIFORM:return "uniform";case ALFVEN:return "alfven";case DIFFUSION:return "diffusion";case BRIO_WU:return "brio_wu";case OBLIQUE_ALFVEN:return "oblique_alfven";case LINEAR_WAVE:return "linear_wave";case FIELD_LOOP:return "field_loop";case RJ2A:return "rj2a";case ORSZAG_TANG:return "orszag_tang";case BLAST_LOW_BETA:return "blast_low_beta";default:throw std::runtime_error("unknown case id");}}
static void json_array(std::ostream& o,const std::vector<double>& a) {
    o<<'[';for(size_t i=0;i<a.size();++i){if(i)o<<',';o<<std::setprecision(17)<<a[i];}o<<']';
}
class Output {
    const Setup& setup;const World& w;const Mesh& mesh;const HostMesh& host;
    std::ofstream diagnostics;int sequence=0;
    std::ofstream midplane_metrics;
    Buffer<State> sensor_state;Buffer<double> sensor_potential;
    std::unique_ptr<TransportOutput> transport_output;
    bool have_budget_origin=false;
    double initial_mass=0,initial_energy=0,initial_angular=0,initial_mass_out=0,initial_energy_out=0,initial_angular_out=0;
    double initial_floor_mass=0,initial_correction=0,initial_cooling=0,initial_floor_potential=0;
public:
    Output(const Setup& ss,const World& ww,const Mesh& mm,const HostMesh& hh):setup(ss),w(ww),mesh(mm),host(hh) {
        if(w.rank==0) {
            fs::create_directories(setup.output);
            if(fs::exists(fs::path(setup.output)/"metadata.json"))throw std::runtime_error("output already contains metadata.json; select a fresh --output directory");
        }
        if(setup.p.transport_diagnostics)transport_output=std::make_unique<TransportOutput>(setup.output,w,host.view,setup.p);
        if(w.rank)return;
        if(setup.midplane_diagnostics) {
            midplane_metrics.open(fs::path(setup.output)/"midplane_metrics.csv");
            midplane_metrics<<"step,time,mu_r1p1_r1p5,Btheta_squared_mean,raw_pressure_mean,radial_width,delta_psi_r1p1_r1p5\n";
        }
        const auto& p=setup.p;std::ofstream meta(fs::path(setup.output)/"metadata.json");meta<<std::setprecision(17);
        meta<<"{\n\"format_version\":2,\n\"geometry\":\""<<(p.spherical?"spherical":"cartesian")<<"\",\n\"case\":\""<<case_name(p.problem)<<"\",\n\"nx1\":"<<p.nx<<",\"nx2\":"<<p.ny<<",\"gamma\":"<<p.gamma<<",\n";
        meta<<"\"checkpoint_version\":"<<pymhd_checkpoint::version<<",\"restarted\":"<<(setup.restart.empty()?"false":"true")
            <<",\"checkpoint_interval\":"<<setup.checkpoint_interval<<",\"retry_streak_limit\":"<<setup.retry_streak_limit
            <<",\"retry_streak_threshold\":"<<setup.retry_streak_threshold<<",\n";
        meta<<"\"source_fingerprint\":\""<<PYMHD_SOURCE_FINGERPRINT<<"\",\n";
        meta<<"\"wave_mode\":"<<p.wave_mode<<",\"noise_seed\":"<<p.noise_seed<<",\"noise_amplitude\":"<<p.noise_amplitude<<",\n";
        meta<<"\"midplane_diagnostics\":"<<setup.midplane_diagnostics<<",\"midplane_pressure\":\"raw total-energy pressure\",\n";
        meta<<"\"eta\":"<<p.eta<<",\"eta_model\":"<<p.eta_model<<",\"wave_amplitude\":"<<p.wave_amplitude<<",\"cfl\":"<<p.cfl<<",\n";
        meta<<"\"entropy_switch\":"<<p.entropy_switch<<",\"entropy_ratio\":"<<p.entropy_ratio<<",\"entropy_shock_threshold\":"<<p.entropy_shock_threshold<<",\"ct_energy_correction\":"<<p.ct_energy_correction<<",\"max_step_retries\":"<<p.max_step_retries<<",\"reconstruction_cache\":"<<p.reconstruction_cache<<",\n";
        meta<<"\"entropy_contact_threshold\":"<<p.entropy_contact_threshold<<",\n";
        meta<<"\"ohmic_cooling_discretization\":\"signed compatible magnetic-energy exchange; distinct from physical J.eta.J\",\n";
        meta<<"\"metric_ppm\":"<<p.metric_ppm<<",\"diffusivity_normalization\":\"etaP=alpha_m*VA_mid*(2*epsilon*R)*F; strong model uses isothermal sqrt(p/rho)\",\n";
        meta<<"\"retry_flat_same_dt\":"<<p.retry_flat_same_dt<<",\n";
        meta<<"\"negative_pressure_policy\":\""<<(p.negative_pressure_floor?"floor":"retry")<<"\",\n";
        meta<<"\"internal_energy_sync\":\""<<(p.internal_energy_sync==1?"modern":"legacy")<<"\",\n";
        meta<<"\"dual_energy_eta\":"<<p.dual_energy_eta<<",\n";
        meta<<"\"metric_cache\":"<<(mesh.view.radial_geometry?"true":"false")<<",\n";
        meta<<"\"ct_emf_scheme\":\""<<(p.ct_emf_scheme==1?"face_hll":"corner_hll")<<"\",\n";
        meta<<"\"cfl_staggered\":"<<p.cfl_staggered<<",\n";
        meta<<"\"transport_diagnostics\":"<<p.transport_diagnostics<<",\n";
        meta<<"\"inner_normal_field\":\""<<(p.inner_flux_project?"project":"evolve")<<"\",\n";
        meta<<"\"cfl_definition\":\"cfl/max_cell(sum_directional_wave_rates+2*eta_max*sum_inverse_length_squared); not PLUTO's directional-average CFL convention\",\n";
        meta<<"\"ranks\":"<<w.size<<",\"backend\":\"";
#ifdef PYMHD_CUDA
        meta<<"cuda";
#else
        meta<<"cpu";
#endif
        meta<<"\",\"reconstruction\":\""<<(p.metric_ppm?"metric-weighted physical-coordinate PPM":"PPM in volume coordinates")<<"\",\"riemann\":\"HLL\",\"time_integrator\":\"SSP-RK3\",\"divergence_control\":\"UCT-HLL axisymmetric flux-function CT\",\n";
        meta<<"\"initial_equilibrium\":\""<<(p.problem==DISK?(p.initial_equilibrium?"Zanni selfsimilar magnetohydrostatic ODE; SF2014 original initializer ambiguous":"hydrostatic polytrope with non-force-free CT magnetic field; SF2014 prose alternative"):"analytic test initial condition")<<"\",\n";
        meta<<"\"diagnostic_normalization\":\""<<(p.spherical?"upper-domain per radian azimuth; mass flux estimated from adjacent cell":"Cartesian per unit third length")<<"\",\n";
        meta<<"\"x1_min\":"<<p.xlo<<",\"x1_max\":"<<p.xhi<<",\"x2_min\":"<<p.ylo<<",\"x2_max\":"<<p.yhi<<",\n\"x1_faces\":";json_array(meta,mesh.global_xfaces);meta<<",\n\"x2_faces\":";json_array(meta,mesh.global_yfaces);
        meta<<",\n\"gm\":"<<p.gm<<",\"epsilon\":"<<p.h<<",\"mu0\":"<<p.mu0<<",\"mu_ref\":"<<p.mu_ref<<",\"b0\":"<<p.b0<<",\"alpha_m\":"<<p.alpha_m<<",\"chi\":"<<p.chi<<",\"corona_delta\":"<<p.corona<<",\"cool_ohmic\":"<<p.cooling<<"\n}\n";
        diagnostics.open(fs::path(setup.output)/"diagnostics.csv");
        diagnostics<<"step,time,dt,mass,energy,min_rho,min_pressure,max_divB,l1_divB,floor_events,wall_seconds,angular_momentum,mdot_inner,mdot_outer,floor_mass,floor_energy,ct_boundary_energy,entropy_events,entropy_energy,ct_energy,compatible_cooling,physical_ohmic,negative_exchange,internal_events,internal_energy,rejected_steps,flat_steps,cfl_dt,requested_dt,last_step_retries,last_trial_stages,cfl_adv_rate,cfl_diff_rate,cfl_cell_i,cfl_cell_j,cfl_rank,boundary_mass_out,boundary_energy_out,boundary_angular_momentum_out,floor_potential_energy,mass_balance_error,energy_balance_error,angular_momentum_balance_error,dual_unreliable_shock_cells\n";
    }
    void write(Solver& s,int step,double time,double dt,double wall,bool snapshot,double cfl_dt=0,double requested_dt=0) {
        // Report the transport coefficient of the accepted state, rather than
        // a stale table from the preceding RHS stage. This does not change U.
        s.refresh_midplane(s.u.p,s.a.p);
        const auto p=s.p;const auto m=host.view;auto states=s.u.download();auto aa=s.a.download();
        std::vector<State> sensor_states;
        if((p.entropy_switch==3&&p.internal_energy_sync==1)||p.transport_diagnostics) {
            // Accepted U has not exchanged its final-stage ghost cells. Refresh
            // a diagnostic COPY, never project boundaries onto the live state
            // outside the RK accounting. This is a sampled risk indicator,
            // not a count of actual candidate-stage synchronization decisions.
            if(sensor_state.n!=s.u.n){sensor_state.allocate(s.u.n);sensor_potential.allocate(s.a.n);}
            sensor_state.copy_from(s.u);sensor_potential.copy_from(s.a);
            s.exchange(sensor_state.p,sensor_potential.p);sensor_states=sensor_state.download();
        }
        auto mid=p.eta_model?s.mid_global.download():std::vector<Midplane>{};Parameters output_params=p;
        output_params.midplane=mid.empty()?nullptr:mid.data();output_params.primitive_cache=nullptr;
        if(transport_output)transport_output->write(s,states,aa,sensor_states,step,time,snapshot);
        if(setup.midplane_diagnostics)write_midplane_sample(midplane_metrics,w,states,aa,m,p,step,time);
        const int nf=18;std::vector<double> local(static_cast<size_t>(m.nx)*m.ny*nf);double sums[8]={0,0,1e300,1e300,0,0,0,0};double bad=0,unreliable_shocks=0;
        for(int j=0;j<m.ny;++j)for(int i=0;i<m.nx;++i) {
            int ii=i+NG,jj=j+NG;State u=states[m.cell(ii,jj)];Primitive v=primitive(u,m.radius(ii,jj,p),p);
            for(int c=0;c<NV;++c)if(!std::isfinite(u[c]))bad=1;
            double pressure=raw_pressure(u,m.radius(ii,jj,p),p);
            if(!(u[RHO]>0)||!(pressure>0)||!std::isfinite(pressure))bad=1;
            if(p.entropy_switch==3&&p.internal_energy_sync==1&&!modern_energy_reliable(u,m.radius(ii,jj,p),p)
                &&shock_sensor(sensor_states.data(),ii,jj,m,output_params))unreliable_shocks+=1;
            double vol=m.volume(ii,jj,p),div=(m.area1(ii+1,jj,p)*field1(aa.data(),ii+1,jj,m,p)-m.area1(ii,jj,p)*field1(aa.data(),ii,jj,m,p)
                +m.area2(ii,jj+1,p)*field2(aa.data(),ii,jj+1,m,p)-m.area2(ii,jj,p)*field2(aa.data(),ii,jj,m,p))/vol;
            double* row=local.data()+(static_cast<size_t>(j)*m.nx+i)*nf;
            row[0]=m.xc[ii];row[1]=m.yc[jj];row[2]=v[RHO];row[3]=v[M1];row[4]=v[M2];row[5]=v[M3];row[6]=v[ENE];row[7]=v[B1];row[8]=v[B2];row[9]=v[B3];
            row[10]=.25*(aa[m.vertex(ii,jj)]+aa[m.vertex(ii+1,jj)]+aa[m.vertex(ii,jj+1)]+aa[m.vertex(ii+1,jj+1)]);
            row[11]=div;row[12]=vol;row[13]=u[ENTR]/u[RHO];row[14]=u[ENE];row[15]=u[EINT];
            row[16]=eta_at(states.data(),ii,m.xc[ii],m.yc[jj],m,output_params);
            row[17]=(sq(v[B1])+sq(v[B2]))/(2*v[ENE]);
            sums[0]+=u[RHO]*vol;sums[1]+=(u[ENE]+u[RHO]*potential(m.xc[ii],p))*vol;
            sums[2]=std::min(sums[2],v[RHO]);sums[3]=std::min(sums[3],v[ENE]);sums[4]=std::max(sums[4],std::abs(div));sums[5]+=std::abs(div)*vol;sums[6]+=vol;sums[7]+=u[M3]*vol;
        }
        double floors=0;for(double n:s.floors.download())floors+=n;
        if(w.maximum(bad)>0)throw std::runtime_error("invalid raw conserved state at output/final time");
        double repairs[3]={0,0,0};for(double n:s.floor_mass.download())repairs[0]+=n;for(double n:s.floor_energy.download())repairs[1]+=n;for(double n:s.boundary_energy.download())repairs[2]+=n;
        double thermal[8]{};
        for(auto q:s.thermal.download()){thermal[0]+=q.entropy_events;thermal[1]+=q.entropy_energy;thermal[2]+=q.ct_energy;thermal[3]+=q.compatible_cooling;thermal[4]+=q.physical_ohmic;thermal[5]+=q.negative_exchange;thermal[6]+=q.internal_events;thermal[7]+=q.internal_energy;}
        double flux_integrals[4]{};
        for(auto q:s.flux_accounting.download()){flux_integrals[0]+=q.mass_out;flux_integrals[1]+=q.energy_out;flux_integrals[2]+=q.angular_momentum_out;flux_integrals[3]+=q.floor_potential;}
        double mdot[2]={0,0};
        for(int j=NG;j<NG+m.ny;++j) {
            if(w.rank==0){auto v=primitive(states[m.cell(NG,j)],m.radius(NG,j,p),p);mdot[0]-=m.area1(NG,j,p)*v[RHO]*v[M1];}
            if(w.rank==w.size-1){int i=NG+m.nx-1;auto v=primitive(states[m.cell(i,j)],m.radius(i,j,p),p);mdot[1]+=m.area1(i+1,j,p)*v[RHO]*v[M1];}
        }
        std::vector<double> all;std::vector<int>counts(w.size),displs(w.size);int owned=static_cast<int>(local.size());
#ifdef PYMHD_MPI
        MPI_Gather(&owned,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
        if(w.rank==0){for(int r=1;r<w.size;++r)displs[r]=displs[r-1]+counts[r-1];all.resize(displs.back()+counts.back());}
        if(snapshot)MPI_Gatherv(local.data(),owned,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        double combined[8]{};for(int k=0;k<8;++k){MPI_Op op=(k==2||k==3)?MPI_MIN:k==4?MPI_MAX:MPI_SUM;MPI_Reduce(sums+k,combined+k,1,MPI_DOUBLE,op,0,MPI_COMM_WORLD);}std::copy(combined,combined+8,sums);
        double f=0;MPI_Reduce(&floors,&f,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);floors=f;
        double repaired[3]{},fluxes[2]{};MPI_Reduce(repairs,repaired,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);MPI_Reduce(mdot,fluxes,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        std::copy(repaired,repaired+3,repairs);std::copy(fluxes,fluxes+2,mdot);
        double combined_thermal[8]{};MPI_Reduce(thermal,combined_thermal,8,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);std::copy(combined_thermal,combined_thermal+8,thermal);
        double combined_flux[4]{};MPI_Reduce(flux_integrals,combined_flux,4,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);std::copy(combined_flux,combined_flux+4,flux_integrals);
        double combined_shocks=0;MPI_Reduce(&unreliable_shocks,&combined_shocks,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);unreliable_shocks=combined_shocks;
#else
        if(snapshot)all=local;
        counts[0]=owned;
#endif
        if(w.rank)return;
        const double correction=repairs[1]+repairs[2]+thermal[1]+thermal[2]+thermal[7];
        // On restart these are segment balances relative to the loaded state;
        // the cumulative flux/correction columns retain the original history.
        if(!have_budget_origin){have_budget_origin=true;initial_mass=sums[0];initial_energy=sums[1];initial_angular=sums[7];
            initial_mass_out=flux_integrals[0];initial_energy_out=flux_integrals[1];initial_angular_out=flux_integrals[2];
            initial_floor_mass=repairs[0];initial_correction=correction;initial_cooling=thermal[3];initial_floor_potential=flux_integrals[3];}
        const double mass_error=(sums[0]-initial_mass)+(flux_integrals[0]-initial_mass_out)-(repairs[0]-initial_floor_mass);
        const double energy_error=(sums[1]-initial_energy)+(flux_integrals[1]-initial_energy_out)+(thermal[3]-initial_cooling)
            -(correction-initial_correction)-(flux_integrals[3]-initial_floor_potential);
        const double angular_error=(sums[7]-initial_angular)+(flux_integrals[2]-initial_angular_out);
        diagnostics<<std::setprecision(17)<<step<<','<<time<<','<<dt<<','<<sums[0]<<','<<sums[1]<<','<<sums[2]<<','<<sums[3]<<','<<sums[4]<<','<<sums[5]/sums[6]<<','<<floors<<','<<wall<<','<<sums[7]<<','<<mdot[0]<<','<<mdot[1]<<','<<repairs[0]<<','<<repairs[1]<<','<<repairs[2];
        for(double q:thermal)diagnostics<<','<<q;
        diagnostics<<','<<s.rejected_steps<<','<<s.flat_steps<<','<<cfl_dt<<','<<requested_dt<<','<<s.last_step_retries<<','<<s.last_trial_stages
                   <<','<<s.last_rate.adv<<','<<s.last_rate.diff<<','<<(s.last_rate.cell>=0?s.last_rate.cell%p.nx:-1)
                   <<','<<(s.last_rate.cell>=0?s.last_rate.cell/p.nx:-1)<<','<<s.last_rate_rank;
        for(double value:flux_integrals)diagnostics<<','<<value;
        diagnostics<<','<<mass_error<<','<<energy_error<<','<<angular_error<<','<<unreliable_shocks<<'\n';diagnostics.flush();
        std::cout<<"step="<<step<<" t="<<std::setprecision(9)<<time<<" dt="<<dt<<" min(rho,p)="<<sums[2]<<','<<sums[3]<<" divB="<<sums[4]<<" floors="<<floors<<" thermal_dE="<<thermal[1]+thermal[7]<<" ct_dE="<<thermal[2]<<" retries="<<s.rejected_steps<<" last_retries="<<s.last_step_retries<<" dt_over_cfl="<<(cfl_dt>0?dt/cfl_dt:0)<<" wall="<<wall<<"s\n"<<std::flush;
        if(snapshot) {
            std::ostringstream name;name<<"snapshot_"<<std::setfill('0')<<std::setw(6)<<sequence++<<".csv";
            fs::path target=fs::path(setup.output)/name.str(),temp=target.string()+".tmp";std::ofstream out(temp);out<<std::setprecision(17);
            out<<"# time="<<time<<" step="<<step<<"\nx1,x2,rho,vr,vtheta,vphi,pressure,Br,Btheta,Bphi,psi,divB,entropy,total_energy,internal_energy,eta_p,mu_poloidal\n";
            // Gathered rank blocks are radial slabs; output is always global j-major.
            for(int j=0;j<p.ny;++j)for(int r=0;r<w.size;++r) {
                int nx=counts[r]/(p.ny*nf);
                for(int i=0;i<nx;++i){const double* row=all.data()+displs[r]+(static_cast<size_t>(j)*nx+i)*nf;for(int c=0;c<12;++c){if(c)out<<',';out<<row[c];}out<<','<<row[13]<<','<<row[14]<<','<<row[15]<<','<<row[16]<<','<<row[17]<<'\n';}
            }
            out.close();fs::rename(temp,target);
        }
    }
};
