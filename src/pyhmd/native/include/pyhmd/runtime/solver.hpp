#pragma once
#include "pyhmd/numerics/cfl.hpp"
#include "pyhmd/boundaries/physical.hpp"
#include "pyhmd/boundaries/gpu_scan.hpp"
#include "pyhmd/distributed/world.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <functional>

class Solver {
public:
    Parameters p;MeshView m;const World& world;
    Buffer<State> u,base,stage,tmp,rhs,flux1,flux2;
    Buffer<Primitive> primitive_cache,edge_xl,edge_xr,edge_yl,edge_yr;
    Buffer<Current> current_cache;
    Buffer<FaceSignal> signal1,signal2;
    Buffer<EtaGeometry> eta_cell_geometry,eta_radial_geometry,eta_theta_geometry,eta_corner_geometry;
    // floors counts individual density/pressure repair events. floor_mass and
    // floor_energy are cumulative extensive injections per cell, weighted by
    // their final-step SSP RK3 contribution (1/6, 2/3, 1 for the three stages).
    // They already include cell volume; diagnostics should sum, not reweight.
    Buffer<double> a,abase,astage,atmp,arhs,floors,floor_mass,floor_energy,boundary_energy,send_l,send_r,recv_l,recv_r;
    Buffer<Rate> rates,rate_result;
    Buffer<Midplane> mid_local,mid_global;
    Buffer<State> resistive1,resistive2;
    Buffer<double> resistive_a,ohmic_heat,physical_heat,trial_floors,trial_mass,trial_energy,trial_boundary,trial_ct;
    Buffer<ThermalAccounting> thermal,trial_thermal;
    Buffer<FluxAccounting> flux_accounting,trial_flux_accounting;
    Buffer<TransportFace> transport_radial,transport_polar,trial_transport_radial,trial_transport_polar;
    Buffer<TransportVertex> transport_psi,trial_transport_psi;
    Buffer<State> transport_origin_state;
    Buffer<double> transport_origin_psi;
    Buffer<double> periodic_flux_jump;
    bool transport_origin_set=false;
    double transport_origin_time=0;
    long long transport_origin_step=0;
    Buffer<StageFailure> failure;
    Buffer<int> shock_flags;
    Buffer<double> face_velocity1,face_velocity2;
    long long rejected_steps=0,flat_steps=0;
    // Host-side diagnostics: called collectively after a rejected trial and
    // before any accepted accounting is committed. base/abase still contain
    // the last accepted state. A rank with no local failure receives active=0.
    std::function<void(const Solver&,const StageFailure&,int)> rejected_trial_hook;
    StageFailure first_rejected_stage,last_rejected_stage;
    int last_step_retries=0,last_trial_stages=0;
    Rate last_rate;
    int last_rate_rank=0;
    std::string failure_directory;
    int halo_size;
    Solver(Parameters pp,MeshView mm,const World& w):p(pp),m(mm),world(w),u(m.tx*m.ty),base(m.tx*m.ty),stage(m.tx*m.ty),tmp(m.tx*m.ty),rhs(m.tx*m.ty),flux1((m.nx+1)*m.ny),flux2(m.nx*(m.ny+1)),
        primitive_cache(m.tx*m.ty),edge_xl(m.tx*m.ty),edge_xr(m.tx*m.ty),edge_yl(m.tx*m.ty),edge_yr(m.tx*m.ty),
        a((m.tx+1)*(m.ty+1)),abase(a.n),astage(a.n),atmp(a.n),arhs(a.n),floors(m.nx*m.ny),floor_mass(floors.n),floor_energy(floors.n),boundary_energy(floors.n),
        send_l(NV*NG*m.ty+(NG+1)*(m.ty+1)),send_r(send_l.n),recv_l(send_l.n),recv_r(send_l.n),rates((m.nx*m.ny+127)/128),rate_result(1),mid_local(m.nx),mid_global(p.nx) {
        halo_size=static_cast<int>(send_l.n);floors.upload(std::vector<double>(floors.n,0));
        if(world.size>1&&p.lower1==PERIODIC)periodic_flux_jump.allocate(m.ty+1);
        floor_mass.upload(std::vector<double>(floors.n,0));floor_energy.upload(std::vector<double>(floors.n,0));boundary_energy.upload(std::vector<double>(floors.n,0));
        resistive1.allocate(flux1.n);resistive2.allocate(flux2.n);resistive_a.allocate(a.n);
        ohmic_heat.allocate(floors.n);physical_heat.allocate(floors.n);
        trial_floors.allocate(floors.n);trial_mass.allocate(floors.n);trial_energy.allocate(floors.n);trial_boundary.allocate(floors.n);trial_ct.allocate(floors.n);
        thermal.allocate(floors.n);trial_thermal.allocate(floors.n);failure.allocate(1);
        flux_accounting.allocate(floors.n);trial_flux_accounting.allocate(floors.n);
        shock_flags.allocate(m.tx*m.ty);
        face_velocity1.allocate(flux1.n);face_velocity2.allocate(flux2.n);
        if(p.reconstruction_cache&&(p.eta_model!=0||p.eta>0))current_cache.allocate(u.n);
        thermal.upload(std::vector<ThermalAccounting>(floors.n));
        flux_accounting.upload(std::vector<FluxAccounting>(floors.n));
        if(p.transport_diagnostics) {
            transport_radial.allocate(flux1.n);trial_transport_radial.allocate(flux1.n);
            transport_polar.allocate(flux2.n);trial_transport_polar.allocate(flux2.n);
            transport_psi.allocate(a.n);trial_transport_psi.allocate(a.n);
            transport_origin_state.allocate(u.n);transport_origin_psi.allocate(a.n);
            transport_radial.upload(std::vector<TransportFace>(flux1.n));transport_polar.upload(std::vector<TransportFace>(flux2.n));
            transport_psi.upload(std::vector<TransportVertex>(a.n));
        }
        initialize_eta_geometry();
    }
    void ensure_transport_origin(double time=0,long long step=0) {
        if(!p.transport_diagnostics||transport_origin_set)return;
        if(time!=0||step!=0)throw std::runtime_error("transport diagnostics must start from initial state or a complete diagnostic checkpoint");
        transport_origin_state.copy_from(u);transport_origin_psi.copy_from(a);
        transport_origin_time=time;transport_origin_step=step;transport_origin_set=true;
    }
    // nvcc extended lambdas need an addressable enclosing function. A
    // constructor has no function address, so keep these launches in a named
    // public member called once after the constructor allocates its buffers.
    void initialize_eta_geometry() {
        if(p.eta_model!=0) {
            eta_cell_geometry.allocate(u.n);eta_radial_geometry.allocate(flux1.n);
            eta_theta_geometry.allocate(flux2.n);eta_corner_geometry.allocate((m.nx+1)*(m.ny+1));
            auto mv=m;auto params=p;
            auto ec=eta_cell_geometry.p,er=eta_radial_geometry.p,et=eta_theta_geometry.p,ev=eta_corner_geometry.p;
            parallel_for(static_cast<int>(eta_cell_geometry.n),[=]HD(int k){int i=k%mv.tx,j=k/mv.tx;ec[k]=make_eta_geometry(mv.xc[i],mv.yc[j],params);});
            parallel_for(static_cast<int>(eta_radial_geometry.n),[=]HD(int k){int i=k%(mv.nx+1)+NG,j=k/(mv.nx+1)+NG;er[k]=make_eta_geometry(mv.xf[i],mv.yc[j],params);});
            parallel_for(static_cast<int>(eta_theta_geometry.n),[=]HD(int k){int i=k%mv.nx+NG,j=k/mv.nx+NG;et[k]=make_eta_geometry(mv.xc[i],mv.yf[j],params);});
            parallel_for(static_cast<int>(eta_corner_geometry.n),[=]HD(int k){int i=k%(mv.nx+1)+NG,j=k/(mv.nx+1)+NG;ev[k]=make_eta_geometry(mv.xf[i],mv.yf[j],params);});
            p.eta_cell_geometry=ec;p.eta_radial_geometry=er;p.eta_theta_geometry=et;p.eta_corner_geometry=ev;
        }
    }
    void prepare_reconstruction(const State* state,Parameters& pp) {
        // Cache lifetime is one RHS evaluation only. The persistent s.p stays
        // cache-free so boundary projection, stage updates and CFL cannot read
        // primitive values from an earlier state.
        auto mv=m;auto output=primitive_cache.p;Parameters uncached=pp;
        uncached.primitive_cache=nullptr;uncached.edge_xl=uncached.edge_xr=uncached.edge_yl=uncached.edge_yr=nullptr;
        parallel_for(m.tx*m.ty,[=]HD(int k){int i=k%mv.tx,j=k/mv.tx;output[k]=primitive(state[k],mv.radius(i,j,uncached),uncached);});
        pp.primitive_cache=primitive_cache.p;
        auto xl=edge_xl.p,xr=edge_xr.p,yl=edge_yl.p,yr=edge_yr.p;Parameters cell_cached=pp;
        cell_cached.edge_xl=cell_cached.edge_xr=cell_cached.edge_yl=cell_cached.edge_yr=nullptr;
        const int width=m.tx-4,height=m.ty-4;
#ifdef PYMHD_CUDA
        auto cells=primitive_cache.p;
        parallel_for(2*width*height*NV,[=]HD(int k){
            int variable=k%NV,component=k/NV,dir=component/(width*height),loc=component%(width*height);
            int i=loc%width+2,j=loc/width+2,id=mv.cell(i,j);double left,right;
            reconstructed_component_pair(cells,i,j,dir,variable,mv,cell_cached,left,right);
            if(dir==0){xl[id][variable]=left;xr[id][variable]=right;}
            else{yl[id][variable]=left;yr[id][variable]=right;}
        });
#else
        parallel_for(2*width*height,[=]HD(int k){
            int dir=k/(width*height),loc=k%(width*height),i=loc%width+2,j=loc/width+2,id=mv.cell(i,j);
            Primitive left,right;reconstructed_pair(state,i,j,dir,mv,cell_cached,left,right);
            if(dir==0){xl[id]=left;xr[id]=right;}else{yl[id]=left;yr[id]=right;}
        });
#endif
        pp.edge_xl=edge_xl.p;pp.edge_xr=edge_xr.p;pp.edge_yl=edge_yl.p;pp.edge_yr=edge_yr.p;
    }
    void refresh_midplane(const State* state,const double* potential_a) {
        if(p.eta_model==0)return;
        auto mv=m;auto pp=p;auto out=world.size==1?mid_global.p:mid_local.p;
        parallel_for(m.nx,[=]HD(int k){int i=NG+k,j=NG+mv.ny-1;Primitive v=get_primitive(state,i,j,mv,pp);double bp=field2(potential_a,i,NG+mv.ny,mv,pp);Midplane q;q.va=fabs(bp)/sqrt(v[RHO]);q.cs=sqrt(v[ENE]/v[RHO]);q.mu=bp*bp/(2*v[ENE]);out[k]=q;});
        // A single rank owns the complete lookup table already. Keep it on the
        // device instead of downloading and uploading it every RK evaluation.
        if(world.size==1){p.midplane=mid_global.p;return;}
        auto local=mid_local.download();std::vector<Midplane> global(p.nx);
#ifdef PYMHD_MPI
        std::vector<int> count(world.size),disp(world.size);int nr=p.nx/world.size,rem=p.nx%world.size;
        for(int r=0;r<world.size;++r){count[r]=(nr+(r<rem))*static_cast<int>(sizeof(Midplane));disp[r]=(r*nr+std::min(r,rem))*static_cast<int>(sizeof(Midplane));}
        MPI_Allgatherv(local.data(),static_cast<int>(local.size()*sizeof(Midplane)),MPI_BYTE,global.data(),count.data(),disp.data(),MPI_BYTE,MPI_COMM_WORLD);
#else
        global=local;
#endif
        mid_global.upload(global);p.midplane=mid_global.p;
    }
    void prepare_current_cache(const State* state,const double* potential_a,Parameters& pp) {
        if(current_cache.n==0&&(p.eta_model!=0||p.eta>0))current_cache.allocate(u.n);
        if(current_cache.n==0)return;
        auto mv=m;auto out=current_cache.p;Parameters uncached=pp;uncached.current_cache=nullptr;
        const int width=m.tx-2,height=m.ty-2;
        parallel_for(width*height,[=]HD(int k){int i=k%width+1,j=k/width+1;out[mv.cell(i,j)]=cell_current(state,potential_a,i,j,mv,uncached);});
        pp.current_cache=out;
    }
    void exchange(State* state,double* potential_a,double projection_weight=0,double* projection_energy=nullptr) {
        bool lo=world.rank==0,hi=world.rank==world.size-1;
        auto mv=m;auto pp=p;
        parallel_for((m.tx+1)*NG,[=]HD(int k){theta_boundary(k,state,potential_a,mv,pp);});
#ifdef PYMHD_MPI
        if(world.size>1) {
            // Psi(x+L,y)-Psi(x,y) carries the net transverse flux. Copying it
            // periodically without this offset corrupts nonzero mean By.
            std::vector<double> periodic_jump;
            if(p.lower1==PERIODIC) {
                auto ends=periodic_flux_jump.p;
                parallel_for(m.ty+1,[=]HD(int j){ends[j]=lo?-potential_a[mv.vertex(NG,j)]:(hi?potential_a[mv.vertex(NG+mv.nx,j)]:0.);});
                auto local=periodic_flux_jump.download();periodic_jump.resize(m.ty+1);
                MPI_Allreduce(local.data(),periodic_jump.data(),m.ty+1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            }
            int left=lo?(p.lower1==PERIODIC?world.size-1:MPI_PROC_NULL):world.rank-1;
            int right=hi?(p.upper1==PERIODIC?0:MPI_PROC_NULL):world.rank+1;
            double* sl=send_l.p;double* sr=send_r.p;
            parallel_for(halo_size,[=]HD(int k) {
                int nstate=NV*NG*mv.ty;
                if(k<nstate){int c=k%NV,ii=(k/NV)%NG,j=k/(NV*NG);sl[k]=state[mv.cell(NG+ii,j)][c];sr[k]=state[mv.cell(mv.nx+ii,j)][c];}
                else {int q=k-nstate,ii=q%(NG+1),j=q/(NG+1);sl[k]=potential_a[mv.vertex(NG+1+ii,j)];sr[k]=potential_a[mv.vertex(NG+mv.nx-NG+ii,j)];}
            });
            auto hl=send_l.download(),hr=send_r.download();std::vector<double>rl(halo_size,0),rr(halo_size,0);
            MPI_Sendrecv(hl.data(),halo_size,MPI_DOUBLE,left,20,rr.data(),halo_size,MPI_DOUBLE,right,20,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Sendrecv(hr.data(),halo_size,MPI_DOUBLE,right,21,rl.data(),halo_size,MPI_DOUBLE,left,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(!periodic_jump.empty())for(int j=0;j<=m.ty;++j)if(periodic_jump[j]!=0)for(int ii=0;ii<=NG;++ii) {
                const int k=NV*NG*m.ty+j*(NG+1)+ii;
                if(lo)rl[k]-=periodic_jump[j];if(hi)rr[k]+=periodic_jump[j];
            }
            recv_l.upload(rl);recv_r.upload(rr);const double* dl=recv_l.p;const double* dr=recv_r.p;
            parallel_for(halo_size,[=]HD(int k) {
                int nstate=NV*NG*mv.ty;
                if(k<nstate){int c=k%NV,ii=(k/NV)%NG,j=k/(NV*NG);if(left!=MPI_PROC_NULL)state[mv.cell(ii,j)][c]=dl[k];if(right!=MPI_PROC_NULL)state[mv.cell(NG+mv.nx+ii,j)][c]=dr[k];}
                else {int q=k-nstate,ii=q%(NG+1),j=q/(NG+1);
                    if(left!=MPI_PROC_NULL&&!(lo&&pp.lower1==PERIODIC&&ii==NG))potential_a[mv.vertex(ii,j)]=dl[k];
                    if(right!=MPI_PROC_NULL&&ii<NG)potential_a[mv.vertex(NG+mv.nx+1+ii,j)]=dr[k];}
            });
            // Lower-rank owner supplies the duplicate CT interface vertex.
            std::vector<double> shared_send(m.ty+1),shared_recv(m.ty+1);
            for(int j=0;j<=m.ty;++j)shared_send[j]=hr[NV*NG*m.ty+j*(NG+1)+NG];
            MPI_Sendrecv(shared_send.data(),m.ty+1,MPI_DOUBLE,right,22,shared_recv.data(),m.ty+1,MPI_DOUBLE,left,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(lo&&!periodic_jump.empty())for(int j=0;j<=m.ty;++j)if(periodic_jump[j]!=0)shared_recv[j]-=periodic_jump[j];
            if(left!=MPI_PROC_NULL&&!(lo&&p.lower1==PERIODIC)) {
                std::vector<double> values(halo_size,0);std::copy(shared_recv.begin(),shared_recv.end(),values.begin());recv_l.upload(values);const double* vv=recv_l.p;
                parallel_for(m.ty+1,[=]HD(int j){potential_a[mv.vertex(NG,j)]=vv[j];});
            }
        }
#endif
        bool physical_lo=lo&&(world.size==1||p.lower1!=PERIODIC);
        bool physical_hi=hi&&(world.size==1||p.upper1!=PERIODIC);
        parallel_for(2*NG*m.ty,[=]HD(int k){int i=k%(2*NG),j=k/(2*NG);if(i>=NG)i+=mv.tx-2*NG;radial_boundary(mv.cell(i,j),state,mv,pp,physical_lo,physical_hi);});
        parallel_for(m.ty+1,[=]HD(int j){potential_radial_boundary(j,potential_a,mv,pp,physical_lo,physical_hi);});
        if(lo&&p.lower1==PAPER_INNER) {
#ifdef PYMHD_CUDA
            launch_paper_inner_potential(potential_a,mv,pp.inner_flux_project?NG:NG-1);
#else
            paper_inner_potential(0,potential_a,mv,pp);
#endif
        }
        parallel_for((m.tx+1)*NG,[=]HD(int k){theta_boundary(k,state,potential_a,mv,pp);});
        auto be=projection_energy?projection_energy:boundary_energy.p;
        parallel_for(m.tx*m.ty,[=]HD(int k){update_cell_fields(k,state,potential_a,mv,pp,1,be,projection_weight);});
    }
    double timestep() {
        refresh_midplane(u.p,a.p);
        auto mv=m;auto pp=p;auto uu=u.p;auto aa=a.p;Rate result;
#ifdef PYMHD_CUDA
        rate_reduce<<<static_cast<unsigned>(rates.n),128>>>(uu,aa,rates.p,mv,pp);cuda_check(cudaGetLastError(),"rate reduction");
        rate_reduce_final<<<1,128>>>(rates.p,static_cast<int>(rates.n),rate_result.p);cuda_check(cudaGetLastError(),"final rate reduction");
        result=rate_result.download()[0];
#else
        for(int k=0;k<m.nx*m.ny;++k)result=combine_rates(result,rate_at(k,uu,aa,mv,pp));
#endif
        if(world.maximum(result.bad)>0)throw std::runtime_error("nonfinite state or nonpositive raw density/pressure detected");
#ifdef PYMHD_MPI
        struct RateRank {double rate;int rank;};RateRank local{result.total,world.rank},winning;
        MPI_Allreduce(&local,&winning,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
        last_rate_rank=winning.rank;
        MPI_Bcast(&result,static_cast<int>(sizeof(Rate)),MPI_BYTE,last_rate_rank,MPI_COMM_WORLD);
#else
        last_rate_rank=0;
#endif
        last_rate=result;double rate=result.total;
        if(!std::isfinite(rate)||rate<=0)throw std::runtime_error("invalid CFL rate");
        return p.cfl/rate;
    }
    void evaluate(State* state,double* avec,double projection_weight=0,double* projection_energy=nullptr,double transport_weight=0) {
        exchange(state,avec,projection_weight,projection_energy);refresh_midplane(state,avec);auto mv=m;auto pp=p;
        if(p.reconstruction_cache)prepare_reconstruction(state,pp);
        if(p.reconstruction_cache)prepare_current_cache(state,avec,pp);
        if(p.entropy_switch&&!(p.entropy_switch==3&&p.internal_energy_sync==1)) {
            auto flags=shock_flags.p;int width=m.tx-2;
            parallel_for(width*(m.ty-2),[=]HD(int k){int i=k%width+1,j=k/width+1;flags[mv.cell(i,j)]=entropy_exclusion_sensor(state,i,j,mv,pp);});
        }
        auto f1=flux1.p,f2=flux2.p,rr=rhs.p;auto ra=arhs.p;auto rf1=resistive1.p,rf2=resistive2.p;auto rra=resistive_a.p,heat=ohmic_heat.p,physical=physical_heat.p;
        auto vel1=face_velocity1.p,vel2=face_velocity2.p;
        if(p.ct_emf_scheme==1&&signal1.n==0){signal1.allocate(u.n);signal2.allocate(u.n);}
        auto sx=p.ct_emf_scheme==1?signal1.p:nullptr,sy=p.ct_emf_scheme==1?signal2.p:nullptr;
        auto tr=trial_transport_radial.p,tt=trial_transport_polar.p;auto tv=trial_transport_psi.p;
        if(p.reconstruction_cache) {
            parallel_for((m.nx+1)*m.ny,[=]HD(int k){face_flux_impl<true>(k,0,state,avec,f1,mv,pp,rf1,vel1,sx,tr,transport_weight);});
            parallel_for(m.nx*(m.ny+1),[=]HD(int k){face_flux_impl<true>(k,1,state,avec,f2,mv,pp,rf2,vel2,sy,tt,transport_weight);});
        } else {
            parallel_for((m.nx+1)*m.ny,[=]HD(int k){face_flux_impl<false>(k,0,state,avec,f1,mv,pp,rf1,vel1,sx,tr,transport_weight);});
            parallel_for(m.nx*(m.ny+1),[=]HD(int k){face_flux_impl<false>(k,1,state,avec,f2,mv,pp,rf2,vel2,sy,tt,transport_weight);});
        }
        if(p.ct_emf_scheme==1) {
            parallel_for((m.nx+1)*2*NG,[=]HD(int k){ghost_face_signal(k,0,state,avec,sx,mv,pp);});
            parallel_for(2*NG*(m.ny+1),[=]HD(int k){ghost_face_signal(k,1,state,avec,sy,mv,pp);});
            parallel_for((m.nx+1)*(m.ny+1),[=]HD(int k){face_based_emf(k,state,avec,sx,sy,ra,mv,pp,rra,tv,transport_weight);});
        } else if(p.reconstruction_cache)
            parallel_for((m.nx+1)*(m.ny+1),[=]HD(int k){corner_emf_impl<true>(k,state,avec,ra,mv,pp,rra);});
        else parallel_for((m.nx+1)*(m.ny+1),[=]HD(int k){corner_emf_impl<false>(k,state,avec,ra,mv,pp,rra);});
        parallel_for(m.nx*m.ny,[=]HD(int k){cell_rhs(k,state,avec,f1,f2,rr,mv,pp,rf1,rf2,rra,heat,physical,vel1,vel2);});
    }
    void write_failure(StageFailure f,int attempt) const {
        if(!f.active)return;
        std::ostringstream message;message<<"failed RK stage="<<f.stage<<" rank="<<world.rank<<" cell=("<<f.i+m.start-NG<<','<<f.j-NG<<") dt="<<std::setprecision(17)<<f.dt<<" raw_pressure="<<f.pressure<<" reason="<<thermal_failure_name(f.reason);
        std::cerr<<message.str()<<'\n';
        if(failure_directory.empty())return;
        std::ofstream out(failure_directory+"/failed_stage_rank"+std::to_string(world.rank)+".txt");
        out<<message.str()<<"\nattempt="<<attempt<<"\ncomponents=rho,mr,mtheta,Lphi,E,Br,Btheta,Bphi,S,u_aux\n"<<std::setprecision(17);
        out<<"thermal_exclusion="<<f.thermal_exclusion<<"\npre_ct_pressure="<<f.pre_ct_pressure
           <<"\nct_delta_energy="<<f.ct_delta_energy<<"\nauxiliary_pressure="<<f.auxiliary_pressure<<'\n';
        const State states[3]={f.source,f.rhs,f.failed};const char* names[3]={"source","rhs","failed"};
        for(int n=0;n<3;++n){out<<names[n];for(int c=0;c<NV;++c)out<<' '<<states[n][c];out<<'\n';}
    }
    double step(double dt) {
        if(!(dt>0)||!std::isfinite(dt))throw std::runtime_error("step requires a finite positive dt");
        ensure_transport_origin();
        first_rejected_stage=StageFailure{};last_rejected_stage=StageFailure{};last_step_retries=0;last_trial_stages=0;
        base.copy_from(u);abase.copy_from(a);const int original_flat=p.flat_reconstruction;
        auto mv=m;const State* b=base.p;const double* ba=abase.p;auto rr=rhs.p;auto ra=arhs.p;
        auto ff=trial_floors.p,fm=trial_mass.p,fe=trial_energy.p,be=trial_boundary.p,ce=trial_ct.p;
        auto ta=trial_thermal.p;auto fail=failure.p;auto heat=ohmic_heat.p,physical=physical_heat.p;
        auto fb=trial_flux_accounting.p;auto f1=flux1.p,f2=flux2.p;
        auto tr=trial_transport_radial.p,tt=trial_transport_polar.p;auto tv=trial_transport_psi.p;
        const int nr=static_cast<int>(transport_radial.n),nt=static_cast<int>(transport_polar.n),nv=static_cast<int>(transport_psi.n);
        const int nd=std::max(nv,std::max(nr,nt));
        for(int attempt=0;attempt<=p.max_step_retries;++attempt) {
            if(attempt){u.copy_from(base);a.copy_from(abase);p.flat_reconstruction=1;}
            stage.copy_from(u);astage.copy_from(a);tmp.copy_from(u);atmp.copy_from(a);
            parallel_for(m.nx*m.ny,[=]HD(int k){ff[k]=fm[k]=fe[k]=be[k]=ce[k]=0;ta[k]=ThermalAccounting{};fb[k]=FluxAccounting{};});
            parallel_for(1,[=]HD(int){*fail=StageFailure{};});
            if(p.transport_diagnostics)parallel_for(nd,[=]HD(int k){if(k<nr)tr[k]=TransportFace{};if(k<nt)tt[k]=TransportFace{};if(k<nv)tv[k]=TransportVertex{};});
            Parameters pp=p;State* source=u.p;double* pa=a.p;
            StageFailure f;bool bad=false;last_trial_stages=0;
            for(int s=0;s<3;++s) {
                evaluate(source,pa,s==2?2./3:1./6,be,dt*(s==2?2./3:1./6));State* dest=s==0?stage.p:(s==1?tmp.p:u.p);double* da=s==0?astage.p:(s==1?atmp.p:a.p);
                double w0=s==0?0.:(s==1?.75:1./3),w1=1-w0;
                double final_weight=s==0?1./6:(s==1?2./3:1.),rhs_weight=s==2?2./3:1./6;
                Parameters recovery=pp;recovery.shock_flags=shock_flags.p;
                if(p.reconstruction_cache)recovery.primitive_cache=primitive_cache.p;
                parallel_for(m.nx*m.ny,[=]HD(int k){rk_cell(k,b,source,rr,dest,mv,pp,dt,w0,w1,ff,fm,final_weight,fb);});
                parallel_for((m.nx+1)*(m.ny+1),[=]HD(int k){rk_potential(k,ba,pa,ra,da,mv,dt,w0,w1);});
                parallel_for(m.nx*m.ny,[=]HD(int k){
                    int i=k%mv.nx+NG,j=k/mv.nx+NG,id=mv.cell(i,j);
                    double pre_ct_pressure=raw_pressure(dest[id],mv.radius(i,j,pp),pp),energy_before=dest[id][ENE];
                    update_cell_fields(id,dest,da,mv,pp,pp.ct_energy_correction,ce,final_weight);
                    double ct_delta=dest[id][ENE]-energy_before,auxiliary_pressure=(pp.gamma-1)*dest[id][EINT];int reason=THERMAL_OK;
                    bool valid=recover_thermal(k,dest,source,mv,recovery,ff,fe,ta,final_weight,&reason);
                    double raw=raw_pressure(dest[id],mv.radius(i,j,pp),pp);
                    if(!valid)capture_failure(fail,s+1,i,j,dt,source[id],rr[id],dest[id],raw,reason,
                        pp.entropy_switch==3&&pp.internal_energy_sync==1?-1:excluded_thermal_cell(source,i,j,mv,recovery),pre_ct_pressure,ct_delta,auxiliary_pressure);
                    double weight=dt*rhs_weight*mv.volume(i,j,pp);
                    if(pp.cooling)ta[k].compatible_cooling+=weight*heat[k];
                    ta[k].physical_ohmic+=weight*physical[k];
                    ta[k].negative_exchange+=weight*fmax(0.,-heat[k]);
                    accumulate_boundary_flux(k,f1,f2,mv,pp,dt*rhs_weight,fb);
                });
                ++last_trial_stages;
                // Do not evaluate fluxes from a rejected intermediate state.
                // All ranks make the same decision before the next halo exchange.
                f=failure.download()[0];bad=world.maximum(f.active)>0;
                if(bad)break;
                source=dest;pa=da;
            }
            if(!bad) {
                auto totals=thermal.p;auto of=floors.p,om=floor_mass.p,oe=floor_energy.p,ob=boundary_energy.p;
                auto flux_totals=flux_accounting.p;
                parallel_for(m.nx*m.ny,[=]HD(int k){
                    of[k]+=ff[k];om[k]+=fm[k];oe[k]+=fe[k];ob[k]+=be[k];
                    totals[k].entropy_events+=ta[k].entropy_events;totals[k].entropy_energy+=ta[k].entropy_energy;totals[k].ct_energy+=ce[k];
                    totals[k].compatible_cooling+=ta[k].compatible_cooling;totals[k].physical_ohmic+=ta[k].physical_ohmic;totals[k].negative_exchange+=ta[k].negative_exchange;
                    totals[k].internal_events+=ta[k].internal_events;totals[k].internal_energy+=ta[k].internal_energy;
                    flux_totals[k].mass_out+=fb[k].mass_out;flux_totals[k].energy_out+=fb[k].energy_out;
                    flux_totals[k].angular_momentum_out+=fb[k].angular_momentum_out;flux_totals[k].floor_potential+=fb[k].floor_potential;
                });
                if(p.transport_diagnostics) {
                    auto cr=transport_radial.p,ct=transport_polar.p;auto cv=transport_psi.p;
                    parallel_for(nd,[=]HD(int k){
                        if(k<nr)for(int c=0;c<NTF;++c)cr[k].q[c]+=tr[k].q[c];
                        if(k<nt)for(int c=0;c<NTF;++c)ct[k].q[c]+=tt[k].q[c];
                        if(k<nv)for(int c=0;c<NTV;++c)cv[k].q[c]+=tv[k].q[c];
                    });
                }
                if(attempt)++flat_steps;
                p.flat_reconstruction=original_flat;return dt;
            }
            if(attempt==0)first_rejected_stage=f;
            last_rejected_stage=f;last_step_retries=attempt+1;
            if(rejected_trial_hook) {
                try{rejected_trial_hook(*this,f,attempt);}
                catch(...){u.copy_from(base);a.copy_from(abase);p.flat_reconstruction=original_flat;throw;}
            }
            ++rejected_steps;
            if(attempt==p.max_step_retries||!(dt>2e-16)) {
                write_failure(f,attempt);u.copy_from(base);a.copy_from(abase);p.flat_reconstruction=original_flat;
                throw std::runtime_error("RK positivity recovery exhausted; see failed_stage_rank*.txt");
            }
            // First try the more diffusive spatial update at the SAME CFL
            // timestep. Halving before that trial confounds a reconstruction
            // failure with a timestep restriction and needlessly slows runs.
            // Already-flat requests have no distinct same-dt fallback.
            if(!(p.retry_flat_same_dt&&attempt==0&&!original_flat))dt*=.5;
        }
        throw std::runtime_error("unreachable step retry state");
    }
};
