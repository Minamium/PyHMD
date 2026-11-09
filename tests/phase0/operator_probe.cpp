// Dump a frozen common input, one RHS and a fixed accepted-step trajectory.
#include "pyhmd/io/checkpoint.hpp"
#include "pyhmd/problems/initialize.hpp"
#include <filesystem>
#include <iterator>
#include <chrono>
static std::string config_text(const char* path) {
    std::ifstream file(path);if(!file)throw std::runtime_error("probe config missing");
    return {std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
}
static const char* fp_policy() {
#ifdef PYHMD_STRICT_FP
    return "strict";
#else
    return "default";
#endif
}
static void states(const std::filesystem::path& path,const std::vector<State>& values,MeshView m) {
    std::ofstream out(path,std::ios::binary);
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)out.write(reinterpret_cast<const char*>(values[m.cell(i,j)].q),sizeof(State));
    if(!out)throw std::runtime_error("state dump failed");
}
static void vertices(const std::filesystem::path& path,const std::vector<double>& values,MeshView m) {
    std::ofstream out(path,std::ios::binary);
    for(int j=NG;j<=NG+m.ny;++j)for(int i=NG;i<=NG+m.nx;++i)out.write(reinterpret_cast<const char*>(&values[m.vertex(i,j)]),sizeof(double));
    if(!out)throw std::runtime_error("vertex dump failed");
}
int main(int argc,char**argv) {
    World world(argc,argv);
    try {
        if(argc<3)throw std::runtime_error("operator_probe CONFIG OUTPUT [INPUT|-] [STEPS]");
        char key[]="--config";char* setup_argv[]={argv[0],key,argv[1]};auto setup=parse_setup(3,setup_argv);auto p=setup.p;
        Mesh mesh(p,world);HostMesh host(mesh);Solver s(p,mesh.view,world);auto m=mesh.view;
        int count=argc>4?std::stoi(argv[4]):100;double time=0;long long first_step=0;
        std::filesystem::path output=argv[2];std::filesystem::create_directories(output);
        if(world.rank==0){std::ofstream signature(output/"parameters.signature.txt");signature<<pymhd_checkpoint::signature(p);}
        if(argc>3&&std::string(argv[3])!="-"&&std::filesystem::exists(std::filesystem::path(argv[3])/"common_input.bin")) {
            // Probe-only frozen full-domain input, including physical ghosts.
            // This is not a repartitioned restart: only an initial RHS/step test.
            pymhd_checkpoint::Reader input(std::filesystem::path(argv[3])/"common_input.bin");
            if(input.string()!="PyHMD-PROBE-INITIAL-1"||input.string()!=config_text(argv[1])||input.string()!=fp_policy())throw std::runtime_error("probe config/FP policy differs");
            // A test deliberately compares different algorithm/source generations.
            // Record the producer fingerprint; never use this format for restart.
            const auto producer=input.string();(void)producer;
            if(input.value<int>()!=p.nx||input.value<int>()!=p.ny)throw std::runtime_error("probe input grid differs");
            int tx=p.nx+2*NG,ty=p.ny+2*NG;
            auto global=input.buffer<State>(tx*ty);auto ga=input.buffer<double>((tx+1)*(ty+1));input.finish();
            std::vector<State> local(s.u.n);std::vector<double> la(s.a.n);
            for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i)local[m.cell(i,j)]=global[j*tx+m.start+i];
            for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)la[m.vertex(i,j)]=ga[j*(tx+1)+m.start+i];
            s.u.upload(local);s.a.upload(la);s.ensure_transport_origin();
        }
        else if(argc>3&&std::string(argv[3])!="-"){auto clock=pymhd_checkpoint::read(s,argv[3],p.output_interval);time=clock.time;first_step=clock.step;}
        else {initialize(s,host);s.ensure_transport_origin();pymhd_checkpoint::write(s,output/"input",{0,0,p.output_interval,p.output_interval,0});}
        if(world.size==1) {
            pymhd_checkpoint::Writer common(output/"common_input.bin");common.string("PyHMD-PROBE-INITIAL-1");common.string(config_text(argv[1]));common.string(fp_policy());common.string(PYMHD_SOURCE_FINGERPRINT);
            common.value(p.nx);common.value(p.ny);common.buffer(s.u);common.buffer(s.a);common.finish();
        }
        std::ostringstream rank;rank<<std::setw(6)<<std::setfill('0')<<world.rank;std::string suffix="_rank"+rank.str()+".bin";
        auto initial=s.u.download();auto initial_a=s.a.download();states(output/("input_U"+suffix),initial,m);vertices(output/("input_Psi"+suffix),initial_a,m);
        double initial_dt=s.timestep(),rate=s.last_rate.total;
        s.evaluate(s.u.p,s.a.p);synchronize();
        states(output/("rhs"+suffix),s.rhs.download(),m);vertices(output/("emf"+suffix),s.arhs.download(),m);
        // evaluate may fill ghosts; accepted state is reset to the exact input.
        s.u.upload(initial);s.a.upload(initial_a);double dt=0;
        auto begin=std::chrono::steady_clock::now();
        for(int n=0;n<count;++n){dt=s.step(s.timestep());time+=dt;}
        synchronize();double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        double next=(std::floor(time/p.output_interval)+1)*p.output_interval;
        if(next<=time)next+=p.output_interval;
        pymhd_checkpoint::write(s,output/"final",{time,dt,next,p.output_interval,first_step+count});
        states(output/("final_U"+suffix),s.u.download(),m);vertices(output/("final_Psi"+suffix),s.a.download(),m);
        std::ofstream meta(output/("rank"+rank.str()+".json"));
        meta<<std::setprecision(17)<<"{\"rank\":"<<world.rank<<",\"ranks\":"<<world.size<<",\"nx\":"<<m.nx<<",\"ny\":"<<m.ny
            <<",\"start\":"<<m.start<<",\"global_nx\":"<<m.global_nx<<",\"rate\":"<<rate<<",\"initial_dt\":"<<initial_dt
            <<",\"steps\":"<<count<<",\"time\":"<<time<<",\"integration_wall_seconds\":"<<seconds<<",\"rejected_steps\":"<<s.rejected_steps<<",\"source_fingerprint\":\""<<PYMHD_SOURCE_FINGERPRINT<<"\"}\n";
        if(!meta)throw std::runtime_error("probe metadata failed");return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';world.abort(1);return 1;}
}
