// Build this with NVCC/CPU compiler against the same source: the current and
// reconstruction caches must preserve the uncached FV/UCT operator. This is
// an equivalence/timing test, not a test against a physical reference solution.
#define main pymhd_application_main
#include "../src/main.cpp"
#undef main
int main(int argc,char**argv) {
    World w(argc,argv);
    try{auto p=setup_case_defaults(DISK);p.nx=argc>1?std::stoi(argv[1]):1200;p.ny=argc>2?std::stoi(argv[2]):256;
        Mesh mesh(p,w);HostMesh host(mesh);Solver s(p,mesh.view,w);initialize(s,host);
        auto initial=s.u.download();auto a=s.a.download();std::vector<State>rhs;std::vector<double>emf;
        auto original_mesh=s.m;
        for(int cache:{0,1}) {
#ifdef PYMHD_METRIC_CACHE
            s.m=original_mesh;if(!cache){s.m.radial_geometry=nullptr;s.m.polar_geometry=nullptr;}
#endif
            s.p.reconstruction_cache=cache;s.u.upload(initial);s.a.upload(a);s.evaluate(s.u.p,s.a.p);synchronize();
            auto r=s.rhs.download();auto e=s.arhs.download();
            if(cache==0){rhs=r;emf=e;continue;}
            double max_error=0;
            for(int c=0;c<NV;++c){double scale=1e-30,error=0;
                for(int j=NG;j<NG+s.m.ny;++j)for(int i=NG;i<NG+s.m.nx;++i){int id=s.m.cell(i,j);scale=std::max(scale,std::abs(rhs[id][c]));error=std::max(error,std::abs(r[id][c]-rhs[id][c]));}
                max_error=std::max(max_error,error/scale);}
            double scale=1e-30,error=0;for(int j=NG;j<=NG+s.m.ny;++j)for(int i=NG;i<=NG+s.m.nx;++i){int id=s.m.vertex(i,j);scale=std::max(scale,std::abs(emf[id]));error=std::max(error,std::abs(e[id]-emf[id]));}
            max_error=std::max(max_error,error/scale);
            max_error=w.maximum(max_error);
            if(max_error>2e-12)throw std::runtime_error("Cached FV/CT RHS differs from reference");
            std::cout<<"PASS cached/uncached FV and UCT RHS, max component-normalized error="<<max_error<<'\n';
        }
        // Identical frozen initial state for each timing sample; interleave A/B
        // to avoid attributing changing physical dt or state complexity to speed.
        for(int rep=0;rep<5;++rep)for(int c=0;c<2;++c){int cache=(rep+c)%2;s.p.reconstruction_cache=cache;
#ifdef PYMHD_METRIC_CACHE
            s.m=original_mesh;if(!cache){s.m.radial_geometry=nullptr;s.m.polar_geometry=nullptr;}
#endif
            s.u.upload(initial);s.a.upload(a);
            auto begin=std::chrono::steady_clock::now();for(int n=0;n<20;++n)s.evaluate(s.u.p,s.a.p);synchronize();
            double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            std::cout<<"rhs_timing cache="<<cache<<" repeat="<<rep<<" milliseconds="<<seconds*1000/20<<'\n';}
        return 0;
    }catch(const std::exception&e){std::cerr<<"FAIL RHS equivalence: "<<e.what()<<'\n';w.abort(1);return 1;}
}
