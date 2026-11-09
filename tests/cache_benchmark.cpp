// Optional CPU timing/equivalence benchmark, not a GPU performance claim.
#define main pymhd_application_main
#include "../src/main.cpp"
#undef main
int main(int argc,char** argv) {
    World world(argc,argv);auto p=setup_case_defaults(DISK);
    if(argc>1)p.nx=std::stoi(argv[1]);if(argc>2)p.ny=std::stoi(argv[2]);
    const int pairs=argc>3?std::stoi(argv[3]):3;
    Mesh mesh(p,world);HostMesh host(mesh);Solver solver(p,mesh.view,world);initialize(solver,host);
    auto frozen_u=solver.u.download();auto frozen_a=solver.a.download();auto m=mesh.view;
    std::vector<State> reference;std::vector<double> reference_emf;
    auto measure=[&](int cache,bool record) {
        solver.p.reconstruction_cache=cache;solver.u.upload(frozen_u);solver.a.upload(frozen_a);
        auto start=std::chrono::steady_clock::now();solver.evaluate(solver.u.p,solver.a.p);synchronize();
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        auto rhs=solver.rhs.download();auto emf=solver.arhs.download();
        if(reference.empty()){reference=rhs;reference_emf=emf;}
        for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)for(int k=0;k<NV;++k)
            if(rhs[m.cell(i,j)][k]!=reference[m.cell(i,j)][k])throw std::runtime_error("cache changed RHS at identical physical state");
        for(int j=NG;j<=NG+m.ny;++j)for(int i=NG;i<=NG+m.nx;++i)
            if(emf[m.vertex(i,j)]!=reference_emf[m.vertex(i,j)])throw std::runtime_error("cache changed CT RHS at identical physical state");
        if(record)std::cout<<cache<<','<<seconds<<','<<p.nx*p.ny/seconds/1e6<<'\n';
    };
    measure(0,false);measure(1,false);
    std::cout<<"cache,seconds_per_RHS,million_cells_per_second\n"<<std::setprecision(12);
    for(int rep=0;rep<pairs;++rep){measure(rep%2,true);measure(1-rep%2,true);}
}
