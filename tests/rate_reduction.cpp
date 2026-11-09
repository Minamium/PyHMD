// CPU simulation of the two CUDA reduction schedules, including partial blocks.
// Real CUDA compilation and race checking remain required on the GPU system.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
static Rate simulated_reduction(const std::vector<Rate>& input) {
    int blocks=(static_cast<int>(input.size())+127)/128;std::vector<Rate> partial(blocks);
    for(int block=0;block<blocks;++block) {
        Rate lane[128];
        for(int t=0;t<128;++t){int k=block*128+t;if(k<static_cast<int>(input.size()))lane[t]=input[k];}
        for(int stride=64;stride;stride/=2)for(int t=0;t<stride;++t)lane[t]=combine_rates(lane[t],lane[t+stride]);
        partial[block]=lane[0];
    }
    Rate lane[128];
    for(int t=0;t<128;++t)for(int k=t;k<blocks;k+=128)lane[t]=combine_rates(lane[t],partial[k]);
    for(int stride=64;stride;stride/=2)for(int t=0;t<stride;++t)lane[t]=combine_rates(lane[t],lane[t+stride]);
    return lane[0];
}
int main(int argc,char** argv) {
    World world(argc,argv);long cells=0;
    for(int nx:{65,600}) {
        auto p=setup_case_defaults(UNIFORM);p.nx=nx;p.ny=nx==65?17:128;Mesh mesh(p,world);auto m=mesh.view;
        std::vector<State> u(m.tx*m.ty);std::vector<double> a((m.tx+1)*(m.ty+1),0);
        for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){auto v=initial(m.xc[i],m.yc[j],p);v[M1]=.2*sin(.1*i);v[M2]=.1*cos(.2*j);u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
        auto check=[&](double expected_bad) {
            std::vector<Rate> rates(m.nx*m.ny);Rate serial;
            for(int k=0;k<m.nx*m.ny;++k){rates[k]=rate_at(k,u.data(),a.data(),m,p);if(rates[k].total>serial.total||(rates[k].total==serial.total&&(serial.cell<0||rates[k].cell<serial.cell))){double oldbad=serial.bad;serial=rates[k];serial.bad=std::max(serial.bad,oldbad);}serial.bad=std::max(serial.bad,rates[k].bad);}
            auto reduced=simulated_reduction(rates);cells+=static_cast<long>(rates.size());
            if(reduced.total!=serial.total||reduced.bad!=serial.bad||reduced.bad!=expected_bad||reduced.cell!=serial.cell||reduced.adv!=serial.adv||reduced.diff!=serial.diff)throw std::runtime_error("two-stage rate reduction differs from serial winning cell and its correlated rates");
        };
        check(0);u[m.cell(NG+m.nx-1,NG+m.ny-1)][RHO]=0;check(1);
    }
    // The maximum total need not occur at either independent component maximum.
    // A tie across partial blocks must select a stable, real cell record.
    std::vector<Rate> artificial(263);
    for(size_t k=0;k<artificial.size();++k)artificial[k].cell=static_cast<int>(k);
    artificial[3].adv=100;artificial[3].total=100;
    artificial[130].diff=110;artificial[130].total=110;
    artificial[258].adv=60;artificial[258].diff=70;artificial[258].total=130;
    artificial[260]=artificial[258];artificial[260].cell=260;artificial[260].bad=1;
    auto winner=simulated_reduction(artificial);
    if(winner.cell!=258||winner.total!=130||winner.adv!=60||winner.diff!=70||winner.bad!=1)throw std::runtime_error("Correlated rate components or deterministic tie/bad reduction failed");
    std::cout<<"rate reduction simulated cells="<<cells<<" winning cells, correlated adv/diff, ties, and bad flags identical\n";
}
