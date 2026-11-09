#pragma once
#include "pyhmd/runtime/solver.hpp"
#include <array>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <type_traits>
#ifndef PYMHD_SOURCE_FINGERPRINT
#define PYMHD_SOURCE_FINGERPRINT "unversioned-manual-build"
#endif

// Exact conserved/vertex state and accepted accounting. No pointers, caches or
// raw Parameters objects are serialized. Rank decomposition must match.
namespace pymhd_checkpoint {
namespace fs=std::filesystem;
constexpr std::uint32_t version=1, endian_tag=0x01020304;
constexpr std::array<char,16> magic{{'P','y','M','H','D','-','C','H','E','C','K','P','O','I','N','T'}};
struct Clock {double time=0,dt=0,next_output=0,output_interval=0;std::int64_t step=0;};

inline std::string signature(const Parameters& p) {
    // Update the algorithm tag whenever a numerical change breaks restart
    // equivalence. A checkpoint is not a claim of CPU/GPU bitwise equivalence.
    std::ostringstream o;o<<std::setprecision(17)<<"pymhd-nv10-rk3-metric-uct-revision7\n";
    o<<PYMHD_SOURCE_FINGERPRINT<<'\n';
    for(int v:{NV,NG,p.nx,p.ny,p.spherical,p.log_r,p.problem,p.lower1,p.upper1,p.lower2,p.upper2,
               p.cooling,p.eta_model,p.entropy_switch,p.ct_energy_correction,p.max_step_retries,p.reconstruction_cache,
               p.metric_ppm,p.initial_equilibrium,p.retry_flat_same_dt,p.negative_pressure_floor,p.internal_energy_sync,p.ct_emf_scheme,p.cfl_staggered,p.inner_flux_project,p.transport_diagnostics})o<<v<<' ';
    o<<'\n';
    for(double v:{p.xlo,p.xhi,p.ylo,p.yhi,p.gamma,p.cfl,p.gm,p.eta,p.rho_floor,p.p_floor,p.h,p.corona,
                  p.b0,p.mu0,p.bend,p.alpha_m,p.chi,p.mu_ref,p.wave_amplitude,p.entropy_ratio,
                  p.entropy_shock_threshold,p.entropy_contact_threshold,p.dual_energy_eta})o<<v<<' ';
    return o.str();
}
inline std::string shard_name(int rank) {
    std::ostringstream o;o<<"rank"<<std::setw(6)<<std::setfill('0')<<rank<<".bin";return o.str();
}
inline void fail(const std::string& what){throw std::runtime_error("checkpoint: "+what);}
struct Hash {
    std::uint64_t value=14695981039346656037ULL;
    void add(const void* data,std::size_t bytes){const auto* p=static_cast<const unsigned char*>(data);for(std::size_t n=0;n<bytes;++n){value^=p[n];value*=1099511628211ULL;}}
};
struct Writer {
    std::ofstream stream;Hash hash;
    explicit Writer(const fs::path& path):stream(path,std::ios::binary|std::ios::trunc){if(!stream)fail("cannot create "+path.string());}
    void bytes(const void* data,std::size_t count){stream.write(static_cast<const char*>(data),static_cast<std::streamsize>(count));hash.add(data,count);if(!stream)fail("write failed");}
    template<class T>void value(const T& x){static_assert(std::is_trivially_copyable<T>::value,"binary value");bytes(&x,sizeof(x));}
    void string(const std::string& x){value(static_cast<std::uint64_t>(x.size()));bytes(x.data(),x.size());}
    template<class T>void buffer(const Buffer<T>& b){auto data=b.download();value(static_cast<std::uint64_t>(data.size()));bytes(data.data(),sizeof(T)*data.size());}
    void finish(){const auto checksum=hash.value;stream.write(reinterpret_cast<const char*>(&checksum),sizeof(checksum));stream.flush();if(!stream)fail("flush failed");stream.close();if(!stream)fail("close failed");}
};
struct Reader {
    std::ifstream stream;Hash hash;
    explicit Reader(const fs::path& path):stream(path,std::ios::binary){if(!stream)fail("cannot open "+path.string());}
    void bytes(void* data,std::size_t count){stream.read(static_cast<char*>(data),static_cast<std::streamsize>(count));if(!stream)fail("truncated checkpoint");hash.add(data,count);}
    template<class T>T value(){static_assert(std::is_trivially_copyable<T>::value,"binary value");T x;bytes(&x,sizeof(x));return x;}
    std::string string(){auto n=value<std::uint64_t>();if(n>16384)fail("invalid string size");std::string s(static_cast<std::size_t>(n),'\0');bytes(s.data(),s.size());return s;}
    template<class T>std::vector<T> buffer(std::size_t expected){auto n=value<std::uint64_t>();if(n!=expected)fail("buffer layout differs");std::vector<T> data(expected);bytes(data.data(),sizeof(T)*data.size());return data;}
    void finish(){std::uint64_t sum=0;stream.read(reinterpret_cast<char*>(&sum),sizeof(sum));if(!stream)fail("missing checksum");if(sum!=hash.value)fail("checksum mismatch");if(stream.peek()!=std::char_traits<char>::eof())fail("unexpected trailing data");}
};
inline void check_platform() {
    static_assert(std::numeric_limits<double>::is_iec559&&sizeof(double)==8,"checkpoint requires IEEE binary64");
    static_assert(sizeof(State)==NV*sizeof(double),"State must contain only its components");
    static_assert(sizeof(ThermalAccounting)==8*sizeof(double),"accounting layout changed; update checkpoint format");
    static_assert(sizeof(FluxAccounting)==4*sizeof(double),"flux accounting layout changed; update checkpoint format");
    static_assert(sizeof(TransportFace)==NTF*sizeof(double)&&sizeof(TransportVertex)==NTV*sizeof(double),"transport sidecar layout changed; update its schema");
}
inline void write(const Solver& s,const fs::path& directory,const Clock& clock,bool accepted_base=false,
                  const StageFailure* rejected=nullptr,int attempt=-1) {
    check_platform();std::string error;
    try {
        fs::create_directories(directory);const auto target=directory/shard_name(s.world.rank),temporary=fs::path(target.string()+".tmp");
        // Immutable generations prevent an old manifest from advertising a mix
        // of old/new shards if a subsequent write is interrupted.
        if(fs::exists(directory/"checkpoint.meta")||fs::exists(target))fail("destination already contains a checkpoint; use a fresh directory");
        Writer out(temporary);out.bytes(magic.data(),magic.size());out.value(version);out.value(endian_tag);
        out.value(static_cast<std::uint32_t>(NV));out.value(static_cast<std::uint32_t>(NG));
        for(int n:{s.world.rank,s.world.size,s.m.nx,s.m.ny,s.m.tx,s.m.ty,s.m.start,s.m.global_nx})out.value(static_cast<std::int32_t>(n));
        out.string(signature(s.p));out.value(clock.time);out.value(clock.dt);out.value(clock.next_output);out.value(clock.output_interval);out.value(clock.step);
        out.value(static_cast<std::int64_t>(s.rejected_steps));out.value(static_cast<std::int64_t>(s.flat_steps));
        out.buffer(accepted_base?s.base:s.u);out.buffer(accepted_base?s.abase:s.a);
        out.buffer(s.floors);out.buffer(s.floor_mass);out.buffer(s.floor_energy);out.buffer(s.boundary_energy);out.buffer(s.thermal);out.buffer(s.flux_accounting);
        out.finish();fs::rename(temporary,target);
        if(s.p.transport_diagnostics) {
            if(!s.transport_origin_set)fail("transport origin missing");
            const auto extra=directory/("transport_"+shard_name(s.world.rank)),part=fs::path(extra.string()+".tmp");
            Writer transport(part);transport.string("PyMHD-TRANSPORT-1");transport.string(signature(s.p));
            transport.value(out.hash.value);transport.value(clock.time);transport.value(clock.step);
            for(int n:{s.world.rank,s.world.size,s.m.nx,s.m.ny,s.m.start,s.m.global_nx})transport.value(static_cast<std::int32_t>(n));
            transport.value(s.transport_origin_time);transport.value(static_cast<std::int64_t>(s.transport_origin_step));
            transport.buffer(s.transport_radial);transport.buffer(s.transport_polar);transport.buffer(s.transport_psi);
            transport.buffer(s.transport_origin_state);transport.buffer(s.transport_origin_psi);
            transport.finish();fs::rename(part,extra);
        }
        if(rejected) {
            std::ofstream record(directory/("rejected_stage_rank"+std::to_string(s.world.rank)+".txt"));
            const auto& f=*rejected;record<<std::setprecision(17)<<"active="<<f.active<<"\ntime="<<clock.time<<"\nstep="<<clock.step
                <<"\nattempt="<<attempt<<"\nstage="<<f.stage<<"\nrank="<<s.world.rank<<"\ncell="<<f.i+s.m.start-NG<<','<<f.j-NG
                <<"\ntrial_dt="<<f.dt<<"\nreason="<<thermal_failure_name(f.reason)<<"\nraw_pressure="<<f.pressure
                <<"\npre_ct_pressure="<<f.pre_ct_pressure<<"\nct_delta_energy="<<f.ct_delta_energy
                <<"\nauxiliary_pressure="<<f.auxiliary_pressure<<"\nthermal_exclusion="<<f.thermal_exclusion
                <<"\ncomponents=rho,mr,mtheta,Lphi,E,Br,Btheta,Bphi,S,u_aux\n";
            const State states[3]={f.source,f.rhs,f.failed};const char* names[3]={"source","rhs","failed"};
            for(int n=0;n<3;++n){record<<names[n];for(int c=0;c<NV;++c)record<<' '<<states[n][c];record<<'\n';}
            record.flush();if(!record)fail("cannot write rejected-stage record");
        }
    }catch(const std::exception& e){error=e.what();}
    // A failed shard must not leave a manifest advertising a complete restart.
    if(s.world.maximum(error.empty()?0.:1.)>0)fail(error.empty()?"another rank could not write its shard":error);
    if(s.world.rank==0) {
        try {
            const auto target=directory/"checkpoint.meta",temporary=fs::path(target.string()+".tmp");
            std::ofstream meta(temporary);meta<<"PyMHD_CHECKPOINT "<<version<<' '<<s.world.size<<'\n';
            meta<<std::setprecision(17)<<"time "<<clock.time<<" step "<<clock.step<<"\n";
            meta.flush();if(!meta)fail("cannot write manifest");meta.close();if(!meta)fail("cannot close manifest");fs::rename(temporary,target);
        }catch(const std::exception& e){error=e.what();}
    }
    if(s.world.maximum(error.empty()?0.:1.)>0)fail(error.empty()?"root could not commit checkpoint manifest":error);
}
inline Clock read(Solver& s,const fs::path& directory,double output_interval) {
    check_platform();Clock clock;std::string error;
    std::vector<State> state;std::vector<double> potential,floors,mass,energy,boundary;
    std::vector<ThermalAccounting> thermal;std::vector<FluxAccounting> fluxes;std::int64_t rejected=0,flat=0;
    std::vector<TransportFace> transport_r,transport_t;std::vector<TransportVertex> transport_v;
    std::vector<State> transport_origin_state;std::vector<double> transport_origin_psi;
    double transport_origin_time=0;std::int64_t transport_origin_step=0;
    try {
        std::ifstream meta(directory/"checkpoint.meta");std::string name;std::uint32_t format=0;int ranks=0;
        if(!(meta>>name>>format>>ranks)||name!="PyMHD_CHECKPOINT"||format!=version||ranks!=s.world.size)fail("missing, incomplete or incompatible manifest");
        Reader in(directory/shard_name(s.world.rank));std::array<char,16> tag{};in.bytes(tag.data(),tag.size());
        if(tag!=magic||in.value<std::uint32_t>()!=version||in.value<std::uint32_t>()!=endian_tag)fail("unsupported format or byte order");
        if(in.value<std::uint32_t>()!=NV||in.value<std::uint32_t>()!=NG)fail("variable/halo layout differs");
        for(int n:{s.world.rank,s.world.size,s.m.nx,s.m.ny,s.m.tx,s.m.ty,s.m.start,s.m.global_nx})if(in.value<std::int32_t>()!=n)fail("rank/grid layout differs; repartitioned restart is unsupported");
        if(in.string()!=signature(s.p))fail("physics, numerics or algorithm revision differs");
        clock.time=in.value<double>();clock.dt=in.value<double>();clock.next_output=in.value<double>();clock.output_interval=in.value<double>();clock.step=in.value<std::int64_t>();
        if(!std::isfinite(clock.time)||clock.time<0||!std::isfinite(clock.dt)||clock.dt<0||!std::isfinite(clock.next_output)||clock.next_output<=clock.time||clock.output_interval!=output_interval||clock.step<0||clock.step>std::numeric_limits<int>::max())fail("invalid clock or changed output interval");
        rejected=in.value<std::int64_t>();flat=in.value<std::int64_t>();if(rejected<0||flat<0)fail("invalid trial counters");
        state=in.buffer<State>(s.u.n);potential=in.buffer<double>(s.a.n);
        floors=in.buffer<double>(s.floors.n);mass=in.buffer<double>(s.floor_mass.n);energy=in.buffer<double>(s.floor_energy.n);boundary=in.buffer<double>(s.boundary_energy.n);thermal=in.buffer<ThermalAccounting>(s.thermal.n);fluxes=in.buffer<FluxAccounting>(s.flux_accounting.n);
        in.finish();
        if(s.p.transport_diagnostics) {
            Reader transport(directory/("transport_"+shard_name(s.world.rank)));
            if(transport.string()!="PyMHD-TRANSPORT-1"||transport.string()!=signature(s.p))fail("transport schema/signature differs");
            if(transport.value<std::uint64_t>()!=in.hash.value||transport.value<double>()!=clock.time||transport.value<std::int64_t>()!=clock.step)
                fail("transport sidecar belongs to a different checkpoint generation");
            for(int n:{s.world.rank,s.world.size,s.m.nx,s.m.ny,s.m.start,s.m.global_nx})if(transport.value<std::int32_t>()!=n)fail("transport rank/grid differs");
            transport_origin_time=transport.value<double>();transport_origin_step=transport.value<std::int64_t>();
            if(transport_origin_time!=0||transport_origin_step!=0)fail("transport origin must be the initial state");
            transport_r=transport.buffer<TransportFace>(s.transport_radial.n);transport_t=transport.buffer<TransportFace>(s.transport_polar.n);
            transport_v=transport.buffer<TransportVertex>(s.transport_psi.n);transport_origin_state=transport.buffer<State>(s.u.n);
            transport_origin_psi=transport.buffer<double>(s.a.n);transport.finish();
            for(const auto& q:transport_r)for(double x:q.q)if(!std::isfinite(x))fail("nonfinite transport face ledger");
            for(const auto& q:transport_t)for(double x:q.q)if(!std::isfinite(x))fail("nonfinite transport face ledger");
            for(const auto& q:transport_v)for(double x:q.q)if(!std::isfinite(x))fail("nonfinite transport Psi ledger");
            for(double x:transport_origin_psi)if(!std::isfinite(x))fail("nonfinite transport origin Psi");
            for(const auto& q:transport_origin_state)for(double x:q.q)if(!std::isfinite(x))fail("nonfinite transport origin state");
        }
        // Saved state is validated on the host before committing any solver
        // buffer. MeshView contains device pointers in CUDA builds: radius()
        // needs host copies of xc and either cached sin(theta) or theta.
        auto validation_mesh=s.m;
#ifdef PYMHD_CUDA
        std::vector<double> host_xc,host_yc;
        std::vector<PolarGeometry> host_polar;
        if(s.p.spherical) {
            host_xc.resize(s.m.tx);
            cuda_check(cudaMemcpy(host_xc.data(),s.m.xc,sizeof(double)*host_xc.size(),cudaMemcpyDeviceToHost),"checkpoint radial centers");
            validation_mesh.xc=host_xc.data();
            if(s.m.polar_geometry) {
                host_polar.resize(s.m.ty);
                cuda_check(cudaMemcpy(host_polar.data(),s.m.polar_geometry,sizeof(PolarGeometry)*host_polar.size(),cudaMemcpyDeviceToHost),"checkpoint polar geometry");
                validation_mesh.polar_geometry=host_polar.data();
            } else {
                host_yc.resize(s.m.ty);
                cuda_check(cudaMemcpy(host_yc.data(),s.m.yc,sizeof(double)*host_yc.size(),cudaMemcpyDeviceToHost),"checkpoint polar centers");
                validation_mesh.yc=host_yc.data();
            }
        }
#endif
        for(int j=NG;j<NG+s.m.ny;++j)for(int i=NG;i<NG+s.m.nx;++i) {
            auto q=state[s.m.cell(i,j)];for(int c=0;c<NV;++c)if(!std::isfinite(q[c]))fail("nonfinite accepted state");
            if(!(q[RHO]>0)||!(raw_pressure(q,validation_mesh.radius(i,j,s.p),s.p)>0))fail("nonpositive accepted density/pressure");
        }
        for(double value:potential)if(!std::isfinite(value))fail("nonfinite flux function");
    }catch(const std::exception& e){error=e.what();}
    if(s.world.maximum(error.empty()?0.:1.)>0)fail(error.empty()?"another rank rejected its checkpoint":error);
    // Commit only after all ranks validate, without boundary projection or a
    // primitive round-trip that would modify the saved conserved state.
    s.u.upload(state);s.a.upload(potential);s.floors.upload(floors);s.floor_mass.upload(mass);s.floor_energy.upload(energy);s.boundary_energy.upload(boundary);s.thermal.upload(thermal);s.flux_accounting.upload(fluxes);
    if(s.p.transport_diagnostics) {
        s.transport_radial.upload(transport_r);s.transport_polar.upload(transport_t);s.transport_psi.upload(transport_v);
        s.transport_origin_state.upload(transport_origin_state);s.transport_origin_psi.upload(transport_origin_psi);
        s.transport_origin_time=transport_origin_time;s.transport_origin_step=transport_origin_step;s.transport_origin_set=true;
    }
    s.rejected_steps=rejected;s.flat_steps=flat;s.refresh_midplane(s.u.p,s.a.p);return clock;
}
} // namespace pymhd_checkpoint
