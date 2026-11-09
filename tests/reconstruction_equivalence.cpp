// Optimization regression: compare batched PPM with the unchanged scalar
// implementation across both coordinates, all variables, and irregular states.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>
#include <limits>
#include <random>
int main(int argc,char** argv) {
    World world(argc,argv);std::mt19937_64 generator(14073136);
    std::uniform_real_distribution<double> random(0.,1.);double maximum=0,maximum_scaled=0;long count=0;
    for(int spherical=0;spherical<2;++spherical) {
        Parameters p=setup_case_defaults(spherical?DISK:ALFVEN);p.nx=64;p.ny=16;
        Mesh mesh(p,world);auto m=mesh.view;std::vector<State> states(m.tx*m.ty);
        for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
            Primitive v;v[RHO]=std::exp(-12*random(generator));v[ENE]=std::exp(-16*random(generator));
            for(int k=M1;k<=M3;++k)v[k]=4*random(generator)-2;
            for(int k=B1;k<=B3;++k)v[k]=2*random(generator)-1;
            states[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
        }
        Solver solver(p,m,world);solver.u.upload(states);Parameters cached=p;
        solver.prepare_reconstruction(solver.u.p,cached);
        if(solver.p.primitive_cache || solver.p.edge_xl)throw std::runtime_error("persistent cache pointer escaped evaluation");
        for(int j=2;j<m.ty-2;++j)for(int i=2;i<m.tx-2;++i)for(int dir=0;dir<2;++dir)for(int side=0;side<2;++side) {
            auto actual=reconstructed(states.data(),i,j,dir,side,m,p);
            auto from_cache=reconstructed(solver.u.p,i,j,dir,side,m,cached);
            for(int k=0;k<NV;++k) {
                double expected=ppm_edge(states.data(),i,j,dir,side,k,m,p);
                if(k==RHO)expected=std::max(expected,p.rho_floor);
                if(k==ENE)expected=std::max(expected,p.p_floor);
                double difference=std::abs(actual[k]-expected);
                // The scalar and batched expressions can be optimized differently.
                // Bound their rounding difference using the INPUT stencil scale:
                // signed interpolation weights can cancel large entropy samples.
                // Cache/retry/flux comparisons below still require exact equality.
                double scale=std::max(1.,std::abs(expected));
                for(int offset=-2;offset<=2;++offset)
                    scale=std::max(scale,std::abs(get_primitive(states.data(),i+(dir==0?offset:0),j+(dir==1?offset:0),m,p)[k]));
                double tolerance=64*std::numeric_limits<double>::epsilon()*scale;
                maximum=std::max(maximum,difference);++count;
                maximum_scaled=std::max(maximum_scaled,difference/scale);
                if(!std::isfinite(actual[k]) || !std::isfinite(expected) || difference>tolerance) {
                    std::cerr<<std::setprecision(17)<<"PPM mismatch: spherical="<<spherical
                             <<" cell=("<<i<<","<<j<<") dir="<<dir<<" side="<<side
                             <<" component="<<k<<" actual="<<actual[k]<<" expected="<<expected
                             <<" difference="<<difference<<" stencil_scale="<<scale
                             <<" tolerance="<<tolerance<<'\n';
                    throw std::runtime_error("batched PPM changed a reconstructed state");
                }
                if(from_cache[k]!=actual[k])throw std::runtime_error("cached PPM changed a reconstructed state");
                double scalar_left,scalar_right;
                reconstructed_component_pair(cached.primitive_cache,i,j,dir,k,m,cached,scalar_left,scalar_right);
                if((side==0?scalar_left:scalar_right)!=actual[k])throw std::runtime_error("component-parallel PPM changed a reconstructed state");
            }
        }
        cached=p;cached.flat_reconstruction=1;solver.prepare_reconstruction(solver.u.p,cached);
        for(int j=2;j<m.ty-2;++j)for(int i=2;i<m.tx-2;++i) {
            auto actual=reconstructed(solver.u.p,i,j,0,0,m,cached),expected=get_primitive(states.data(),i,j,m,p);
            for(int k=0;k<NV;++k)if(actual[k]!=expected[k])throw std::runtime_error("flat retry used stale PPM cache");
            for(int k=0;k<NV;++k){double l,r;reconstructed_component_pair(cached.primitive_cache,i,j,0,k,m,cached,l,r);if(l!=expected[k]||r!=expected[k])throw std::runtime_error("component-parallel flat retry changed a state");}
        }
        std::vector<double> potential_a(solver.a.n);
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)
            potential_a[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p);
        solver.a.upload(potential_a);cached=p;solver.prepare_reconstruction(solver.u.p,cached);
        solver.prepare_current_cache(solver.u.p,solver.a.p,cached);
        if(cached.current_cache)for(int j=1;j<m.ty-1;++j)for(int i=1;i<m.tx-1;++i) {
            auto actual=cell_current(solver.u.p,solver.a.p,i,j,m,cached),expected=cell_current(solver.u.p,solver.a.p,i,j,m,p);
            for(int d=0;d<3;++d)if(actual[d]!=expected[d])throw std::runtime_error("current cache changed curl B");
        }
        for(int dir=0;dir<2;++dir) {
            int n=dir==0?(m.nx+1)*m.ny:m.nx*(m.ny+1);std::vector<State> original(n),optimized(n);
            for(int k=0;k<n;++k){face_flux(k,dir,states.data(),potential_a.data(),original.data(),m,p);face_flux(k,dir,solver.u.p,solver.a.p,optimized.data(),m,cached);}
            for(int k=0;k<n;++k)for(int q=0;q<NV;++q)
                if(original[k][q]!=optimized[k][q])throw std::runtime_error("cache changed face flux");
        }
        std::vector<double> original_emf(solver.a.n,0),optimized_emf(solver.a.n,0);
        for(int k=0;k<(m.nx+1)*(m.ny+1);++k){corner_emf(k,states.data(),potential_a.data(),original_emf.data(),m,p);corner_emf(k,solver.u.p,solver.a.p,optimized_emf.data(),m,cached);}
        for(size_t k=0;k<original_emf.size();++k)if(original_emf[k]!=optimized_emf[k])throw std::runtime_error("cache changed corner EMF");
    }
    std::cout<<"PPM comparisons="<<count<<" maximum absolute difference="<<maximum
             <<" maximum stencil-scaled difference="<<maximum_scaled<<'\n';
}
