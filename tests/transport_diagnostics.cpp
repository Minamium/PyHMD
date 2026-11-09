// Independent invariants for passive diagnostics, local conservation and retry.
#include "pymhd/setup.hpp"
#include "pymhd/paper_initial.hpp"
#include "pymhd/checkpoint.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>

static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class T>static void equal_buffer(const Buffer<T>& a,const Buffer<T>& b,const char* what) {
    auto x=a.download(),y=b.download();require(x.size()==y.size(),"buffer size differs");
    require(std::memcmp(x.data(),y.data(),x.size()*sizeof(T))==0,what);
}
static void initialize_disk(Solver& s) {
    auto m=s.m;auto p=s.p;
    auto xf=std::vector<double>(),yf=std::vector<double>(),xc=std::vector<double>(),yc=std::vector<double>();
#ifdef PYMHD_CUDA
    xf.resize(m.tx+1);yf.resize(m.ty+1);xc.resize(m.tx);yc.resize(m.ty);
    cuda_check(cudaMemcpy(xf.data(),m.xf,sizeof(double)*xf.size(),cudaMemcpyDeviceToHost),"test xf");
    cuda_check(cudaMemcpy(yf.data(),m.yf,sizeof(double)*yf.size(),cudaMemcpyDeviceToHost),"test yf");
    cuda_check(cudaMemcpy(xc.data(),m.xc,sizeof(double)*xc.size(),cudaMemcpyDeviceToHost),"test xc");
    cuda_check(cudaMemcpy(yc.data(),m.yc,sizeof(double)*yc.size(),cudaMemcpyDeviceToHost),"test yc");
    m.xf=xf.data();m.yf=yf.data();m.xc=xc.data();m.yc=yc.data();m.radial_geometry=nullptr;m.polar_geometry=nullptr;
#endif
    PaperInitial disk(p);std::vector<State> u(s.u.n);std::vector<double>a(s.a.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        auto v=disk.at(m.xc[i],m.yc[j]);
        if(m.yc[j]<0){v[M2]*=-1;v[M3]*=-1;}
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));
        v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));
        u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
}
static void compare_transport(const Solver& a,const Solver& b) {
    equal_buffer(a.transport_radial,b.transport_radial,"radial cumulative ledger differs");
    equal_buffer(a.transport_polar,b.transport_polar,"polar cumulative ledger differs");
    equal_buffer(a.transport_psi,b.transport_psi,"Psi cumulative ledger differs");
    equal_buffer(a.transport_origin_state,b.transport_origin_state,"transport origin state differs");
    equal_buffer(a.transport_origin_psi,b.transport_origin_psi,"transport origin Psi differs");
}
static void hll_cases() {
    Parameters p;for(double velocity:{-10.,0.,10.})for(int direction:{0,1}) {
        Primitive l,r;l[RHO]=1.02;r[RHO]=1;l[ENE]=.01;r[ENE]=.02;
        l[M1+direction]=velocity;r[M1+direction]=velocity*.9;
        l[M3]=.8;r[M3]=.7;l[B3]=.13;r[B3]=-.07;
        l[B1]=.04;r[B1]=.06;l[B2]=.02;r[B2]=.05;
        HllTransport parts;auto plain=hll(l,r,direction,.03,p);
        auto f=hll(l,r,direction,.03,p,nullptr,nullptr,&parts);
        require(std::memcmp(&plain,&f,sizeof(State))==0,"HLL diagnostic changes numerical flux");
        require(std::abs(f[RHO]-parts.mass_adv-parts.mass_jump)<2e-14*(1+std::abs(f[RHO])),"HLL mass split does not close");
        require(std::abs(f[M3]-parts.angular_adv-parts.angular_mag-parts.angular_jump)<2e-14*(1+std::abs(f[M3])),"HLL torque split does not close");
    }
}
static void conservation(Solver& s,const Mesh& mesh) {
    auto m=s.m;auto xf=mesh.xf.download(),yf=mesh.yf.download(),xc=mesh.xc.download(),yc=mesh.yc.download();
    m.xf=xf.data();m.yf=yf.data();m.xc=xc.data();m.yc=yc.data();m.radial_geometry=nullptr;m.polar_geometry=nullptr;
    auto u=s.u.download(),u0=s.transport_origin_state.download();auto a=s.a.download(),a0=s.transport_origin_psi.download();
    auto fr=s.transport_radial.download(),ft=s.transport_polar.download();auto pv=s.transport_psi.download();
    auto fm=s.floor_mass.download(),fe=s.floor_energy.download(),be=s.boundary_energy.download();auto th=s.thermal.download();auto fp=s.flux_accounting.download();
    for(int patch:{0,1}) {
        double delta[3]{},flow[3]{},scale[3]{},correction=0,mass_floor=0;
        int i0=patch?NG+3:NG,i1=patch?NG+m.nx-2:NG+m.nx,j0=patch?NG+5:NG,j1=patch?NG+m.ny-1:NG+m.ny;
        for(int j=j0;j<j1;++j)for(int i=i0;i<i1;++i) {
            int id=m.cell(i,j),k=(j-NG)*m.nx+i-NG;double vol=m.volume(i,j,s.p),phi=potential(m.xc[i],s.p);
            double now[3]={u[id][RHO],u[id][ENE]+phi*u[id][RHO],u[id][M3]};
            double initial[3]={u0[id][RHO],u0[id][ENE]+phi*u0[id][RHO],u0[id][M3]};
            for(int c=0;c<3;++c){delta[c]+=(now[c]-initial[c])*vol;scale[c]+=std::abs(initial[c])*vol;}
            mass_floor+=fm[k];correction+=fe[k]+be[k]+fp[k].floor_potential+th[k].ct_energy+th[k].entropy_energy+th[k].internal_energy-th[k].compatible_cooling;
            int components[3]={TF_MASS,TF_ENERGY,TF_ANGULAR};
            for(int c=0;c<3;++c) {
                int v=components[c],ii=i-NG,jj=j-NG;
                if(i==i0)flow[c]-=fr[jj*(m.nx+1)+ii].q[v];
                if(i==i1-1)flow[c]+=fr[jj*(m.nx+1)+ii+1].q[v];
                if(j==j0)flow[c]-=ft[jj*m.nx+ii].q[v];
                if(j==j1-1)flow[c]+=ft[(jj+1)*m.nx+ii].q[v];
            }
        }
        require(std::abs(delta[0]+flow[0]-mass_floor)<1e-10*(1+scale[0]),"local mass budget does not close");
        require(std::abs(delta[1]+flow[1]-correction)<1e-10*(1+scale[1]),"local total energy budget does not close");
        require(std::abs(delta[2]+flow[2])<1e-10*(1+scale[2]),"local angular momentum budget does not close");
        if(!patch) {
            double original[3]{};for(auto q:fp){original[0]+=q.mass_out;original[1]+=q.energy_out;original[2]+=q.angular_momentum_out;}
            for(int c=0;c<3;++c)require(std::abs(flow[c]-original[c])<1e-11*(1+std::abs(original[c])),"new boundary ledger differs from original accounting");
        }
    }
    double max_psi=0;
    for(int j=NG;j<=NG+m.ny;++j)for(int i=NG;i<=NG+m.nx;++i) {
        int v=m.vertex(i,j);double sum=0;for(double x:pv[v].q)sum+=x;
        max_psi=std::max(max_psi,std::abs(a[v]-a0[v]-sum));
    }
    require(max_psi<1e-11,"actual Psi change does not equal integrated CT terms");
    for(auto* f:{&fr,&ft})for(auto q:*f) {
        require(std::abs(q.q[TF_MASS]-q.q[TF_MASS_ADV]-q.q[TF_MASS_JUMP])<1e-11*(1+std::abs(q.q[TF_MASS])),"integrated mass decomposition does not close");
        require(std::abs(q.q[TF_ANGULAR]-q.q[TF_ANG_ADV]-q.q[TF_ANG_MAG]-q.q[TF_ANG_JUMP])<1e-11*(1+std::abs(q.q[TF_ANGULAR])),"integrated torque decomposition does not close");
    }
    std::cout<<"PASS local/domain M/E/L and original boundary budgets; Psi residual="<<max_psi<<'\n';
}
static void retry_tests(World& w) {
    auto p=setup_case_defaults(UNIFORM);p.nx=32;p.ny=4;p.spherical=p.log_r=0;p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;
    p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;p.gm=p.eta=0;p.eta_model=p.cooling=p.entropy_switch=p.ct_energy_correction=0;
    p.transport_diagnostics=1;p.ct_emf_scheme=1;
    Mesh mesh(p,w);
    auto initialize=[](Solver& s,bool cold){auto m=s.m;std::vector<State> u(s.u.n);std::vector<double>a(s.a.n,0);
        auto xc=std::vector<double>(m.tx);
#ifdef PYMHD_CUDA
        cuda_check(cudaMemcpy(xc.data(),m.xc,sizeof(double)*xc.size(),cudaMemcpyDeviceToHost),"retry xc");
#else
        std::copy(m.xc,m.xc+m.tx,xc.begin());
#endif
        for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=cold?1e-10:.1;v[M1]=cold?.1*(xc[i]-.5):(xc[i]<.5?-2:2);u[m.cell(i,j)]=conserved(v,1,s.p);}
        s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);};
    for(bool cold:{false,true}) {
        Solver retried(p,mesh.view,w);initialize(retried,cold);double dt=retried.timestep()*(cold?1:8);double accepted=retried.step(dt);
        require(retried.rejected_steps>0,"retry fixture did not reject");
        auto flat=p;flat.flat_reconstruction=1;Solver clean(flat,mesh.view,w);initialize(clean,cold);require(clean.step(accepted)==accepted,"clean retry reference rejected");
        equal_buffer(retried.u,clean.u,"retry differs from clean state");equal_buffer(retried.a,clean.a,"retry differs from clean Psi");compare_transport(retried,clean);
        if(cold)require(accepted==dt,"same-dt flat retry changed dt");else require(accepted<dt,"oversized retry did not halve dt");
    }
    std::cout<<"PASS rejected trial exclusion for same-dt flat fallback and dt halving\n";
}
int main(int argc,char**argv){try {
    World w(argc,argv);require(w.size==1,"transport unit test requires one rank; MPI invariants are tested by CLI preflight");hll_cases();
    auto p=setup_case_defaults(DISK);p.nx=64;p.ny=32;Mesh mesh(p,w);Solver plain(p,mesh.view,w);initialize_disk(plain);
    p.transport_diagnostics=1;Solver full(p,mesh.view,w),resumed(p,mesh.view,w);initialize_disk(full);
    auto directory=std::filesystem::temp_directory_path()/("pymhd_transport_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    double time=0,dt=0;
    for(int n=1;n<=100;++n) {
        dt=full.timestep();require(dt==plain.timestep(),"diagnostics changed CFL");
        require(full.step(dt)==plain.step(dt),"diagnostics changed accepted dt");time+=dt;
        equal_buffer(full.u,plain.u,"diagnostics changed conserved bits");equal_buffer(full.a,plain.a,"diagnostics changed Psi bits");
        if(n==50){pymhd_checkpoint::write(full,directory/"half",{time,dt,100,100,n});pymhd_checkpoint::read(resumed,directory/"half",100);}
        else if(n>50)require(resumed.step(resumed.timestep())==dt,"resumed dt differs");
    }
    compare_transport(full,resumed);equal_buffer(full.u,resumed.u,"restart changed state");equal_buffer(full.a,resumed.a,"restart changed Psi");
    pymhd_checkpoint::write(full,directory/"full",{time,dt,100,100,100});pymhd_checkpoint::write(resumed,directory/"resumed",{time,dt,100,100,100});
    for(const char* name:{"rank000000.bin","transport_rank000000.bin","checkpoint.meta"}) {
        std::ifstream a(directory/"full"/name,std::ios::binary),b(directory/"resumed"/name,std::ios::binary);
        require(std::string(std::istreambuf_iterator<char>(a),{})==std::string(std::istreambuf_iterator<char>(b),{}),"full versus 50+50 checkpoint bytes differ");
    }
    conservation(full,mesh);retry_tests(w);
    auto before=resumed.u.download();auto sidecar=directory/"full"/"transport_rank000000.bin";
    std::filesystem::copy_file(directory/"half"/"transport_rank000000.bin",sidecar,std::filesystem::copy_options::overwrite_existing);
    bool rejected=false;try{pymhd_checkpoint::read(resumed,directory/"full",100);}catch(const std::exception&){rejected=true;}
    require(rejected,"mismatched-generation sidecar accepted");
    std::filesystem::remove(sidecar);rejected=false;try{pymhd_checkpoint::read(resumed,directory/"full",100);}catch(const std::exception&){rejected=true;}
    require(rejected,"missing sidecar accepted");
    std::filesystem::copy_file(directory/"resumed"/"transport_rank000000.bin",sidecar);
    {std::fstream file(sidecar,std::ios::in|std::ios::out|std::ios::binary);file.seekg(-16,std::ios::end);char byte=0;file.read(&byte,1);byte^=0x20;file.seekp(-16,std::ios::end);file.write(&byte,1);}
    rejected=false;try{pymhd_checkpoint::read(resumed,directory/"full",100);}catch(const std::exception&){rejected=true;}
    require(rejected,"corrupt sidecar accepted");auto after=resumed.u.download();require(std::memcmp(before.data(),after.data(),sizeof(State)*before.size())==0,"failed sidecar load mutated live state");
    std::filesystem::remove_all(directory);
    std::cout<<"PASS 100-step diagnostic on/off bit identity; exact 50+50 restart; sidecar missing/corrupt/generation rejection\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL transport diagnostics: "<<e.what()<<'\n';return 1;}}
