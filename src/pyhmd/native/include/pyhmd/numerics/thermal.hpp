#pragma once
#include "pyhmd/numerics/finite_volume.hpp"

struct FluxAccounting {
    double mass_out=0,energy_out=0,angular_momentum_out=0,floor_potential=0;
};
// Extensive, time-integrated numerical fluxes through the physical domain.
// E_total includes rho*Phi(cell); radial fluxes therefore use Phi(face).
HD inline void accumulate_boundary_flux(int k,const State* f1,const State* f2,MeshView m,Parameters p,double weight,FluxAccounting* budget) {
    int ii=k%m.nx,jj=k/m.nx,i=ii+NG,j=jj+NG;FluxAccounting& b=budget[k];
    if(m.start+ii==0){State f=f1[jj*(m.nx+1)];double scale=-weight*m.area1(i,j,p);b.mass_out+=scale*f[RHO];b.energy_out+=scale*(f[ENE]+potential(m.xf[i],p)*f[RHO]);b.angular_momentum_out+=scale*f[M3];}
    if(m.start+ii==m.global_nx-1){State f=f1[jj*(m.nx+1)+m.nx];double scale=weight*m.area1(i+1,j,p);b.mass_out+=scale*f[RHO];b.energy_out+=scale*(f[ENE]+potential(m.xf[i+1],p)*f[RHO]);b.angular_momentum_out+=scale*f[M3];}
    if(jj==0){State f=f2[ii];double scale=-weight*m.area2(i,j,p);b.mass_out+=scale*f[RHO];b.energy_out+=scale*(f[ENE]+potential(m.xc[i],p)*f[RHO]);b.angular_momentum_out+=scale*f[M3];}
    if(jj==m.ny-1){State f=f2[m.ny*m.nx+ii];double scale=weight*m.area2(i,j+1,p);b.mass_out+=scale*f[RHO];b.energy_out+=scale*(f[ENE]+potential(m.xc[i],p)*f[RHO]);b.angular_momentum_out+=scale*f[M3];}
}
HD inline void rk_cell(int k,const State* base,const State* in,const State* rhs,State* out,MeshView m,Parameters p,double dt,double w0,double w1,double* floors,double* floor_mass,double final_weight,FluxAccounting* flux_accounting=nullptr) {
    int i=k%m.nx+NG,j=k/m.nx+NG,id=m.cell(i,j);State q;
    for(int c=0;c<NV;++c)q[c]=w0*base[id][c]+w1*(in[id][c]+dt*rhs[id][c]);
    // Pressure is repaired only after CT installs the authoritative B. Never
    // turn nonfinite density into a seemingly valid state.
    if(finite_value(q[RHO]) && q[RHO]>0 && q[RHO]<p.rho_floor) {
        double added_mass=final_weight*(p.rho_floor-q[RHO])*m.volume(i,j,p);
        floor_mass[k]+=added_mass;
        if(flux_accounting)flux_accounting[k].floor_potential+=added_mass*potential(m.xc[i],p);
        q[RHO]=p.rho_floor;floors[k]+=1;
    }
    out[id]=q;
}
HD inline double raw_pressure(State q,double radius,Parameters p) {
    if(!(q[RHO]>0) || !finite_value(q[RHO]) || !(fabs(radius)>0))return NAN;
    const double kinetic=.5*(sq(q[M1])+sq(q[M2])+sq(q[M3]/radius))/q[RHO];
    return (p.gamma-1)*(q[ENE]-kinetic-.5*(sq(q[B1])+sq(q[B2])+sq(q[B3])));
}
HD inline bool modern_energy_reliable(State q,double radius,Parameters p) {
    if(!(q[RHO]>0))return false;
    double internal=raw_pressure(q,radius,p)/(p.gamma-1);
    double kinetic2=(sq(q[M1])+sq(q[M2])+sq(q[M3]/radius))/q[RHO];
    double magnetic2=sq(q[B1])+sq(q[B2])+sq(q[B3]);
    double sound2_density=p.gamma*(p.gamma-1)*internal;
    return sound2_density>p.dual_energy_eta*kinetic2
        &&sound2_density>p.dual_energy_eta*magnetic2&&internal>.5*q[EINT];
}
// Construct a representable energy, then check the same expression used by
// the health validator. One nextafter after a cancelled addition is not enough.
HD inline bool set_pressure(State& q,double radius,double pressure,Parameters p) {
    if(!(q[RHO]>0)||!(pressure>0)||!finite_value(pressure))return false;
    double kinetic=.5*(sq(q[M1])+sq(q[M2])+sq(q[M3]/radius))/q[RHO];
    double magnetic=.5*(sq(q[B1])+sq(q[B2])+sq(q[B3]));
    q[ENE]=kinetic+magnetic+pressure/(p.gamma-1);
    for(int n=0;n<16&&finite_value(q[ENE]);++n) {
        double actual=raw_pressure(q,radius,p);
        if(finite_value(actual)&&actual>=pressure*(1-1e-12))return true;
        q[ENE]=::nextafter(q[ENE],static_cast<double>(INFINITY));
    }
    return false;
}
HD inline void post_ct_pressure_floor(int k,State* u,MeshView m,Parameters p,double* floors,double* floor_energy,double final_weight) {
    int i=k%m.nx+NG,j=k/m.nx+NG,id=m.cell(i,j);State q=u[id];
    for(int c=0;c<NV;++c)if(!finite_value(q[c]))return;
    const double pressure=raw_pressure(q,m.radius(i,j,p),p);
    if(!finite_value(pressure))return;
    if(pressure<p.p_floor) {
        const double old_energy=q[ENE];
        if(set_pressure(q,m.radius(i,j,p),p.p_floor,p)) {
            floor_energy[k]+=final_weight*(q[ENE]-old_energy)*m.volume(i,j,p);
            floors[k]+=1;
        }
        u[id]=q;
    }
}
struct ThermalAccounting {
    double entropy_events=0,entropy_energy=0,ct_energy=0;
    double compatible_cooling=0,physical_ohmic=0,negative_exchange=0;
    double internal_events=0,internal_energy=0;
};
HD inline bool shock_sensor(const State* source,int i,int j,MeshView m,Parameters p) {
    Primitive c=get_primitive(source,i,j,m,p),l=get_primitive(source,i-1,j,m,p),r=get_primitive(source,i+1,j,m,p),d=get_primitive(source,i,j-1,m,p),u=get_primitive(source,i,j+1,m,p);
    double divv=(r[M1]-l[M1])/(m.xc[i+1]-m.xc[i-1])+(u[M2]-d[M2])/((m.yc[j+1]-m.yc[j-1])*(p.spherical?m.xc[i]:1.));
    if(p.spherical)divv+=(2*c[M1]+m.cosine_center(j)/m.sine_center(j)*c[M2])/m.xc[i];
    double low=fmin(c[ENE],fmin(fmin(l[ENE],r[ENE]),fmin(d[ENE],u[ENE])));
    double jump=fmax(fmax(fabs(r[ENE]-c[ENE]),fabs(l[ENE]-c[ENE])),fmax(fabs(u[ENE]-c[ENE]),fabs(d[ENE]-c[ENE])))/fmax(low,p.p_floor);
    return divv<0&&jump>p.entropy_shock_threshold;
}
HD inline bool contact_sensor(const State* source,int i,int j,MeshView m,Parameters p) {
    Primitive c=get_primitive(source,i,j,m,p);
    const int di[4]={-1,1,0,0},dj[4]={0,0,-1,1};
    for(int n=0;n<4;++n) {
        Primitive q=get_primitive(source,i+di[n],j+dj[n],m,p);
        double drho=fabs(c[RHO]-q[RHO])/fmax(c[RHO],q[RHO]);
        double dp=fabs(c[ENE]-q[ENE])/fmax(c[ENE],q[ENE]);
        double db2=0,bscale2=0,dv2=0,vscale2=0;
        for(int a=0;a<3;++a) {
            db2+=sq(c[B1+a]-q[B1+a]);bscale2+=sq(c[B1+a])+sq(q[B1+a]);
            dv2+=sq(c[M1+a]-q[M1+a]);vscale2+=sq(c[M1+a])+sq(q[M1+a]);
        }
        vscale2+=p.gamma*(c[ENE]/c[RHO]+q[ENE]/q[RHO])
            +(sq(c[B1])+sq(c[B2])+sq(c[B3]))/c[RHO]+(sq(q[B1])+sq(q[B2])+sq(q[B3]))/q[RHO];
        double db=sqrt(db2/fmax(bscale2,1e-60));
        double dv=sqrt(dv2/fmax(vscale2,1e-60));
        // A pressure-balanced density contact mixes rho and S nonlinearly.
        // Keep total energy authoritative, even in a magnetically dominated
        // cell. This is a contact safeguard, not a shock-heating criterion.
        if(drho>1e-6&&dp<p.entropy_shock_threshold&&dp<p.entropy_contact_threshold*drho
           &&db<p.entropy_contact_threshold*drho&&dv<p.entropy_contact_threshold*drho)return true;
    }
    return false;
}
HD inline int entropy_exclusion_sensor(const State* source,int i,int j,MeshView m,Parameters p) {
    return (shock_sensor(source,i,j,m,p)?1:0)|((p.entropy_switch!=3&&contact_sensor(source,i,j,m,p))?2:0);
}
HD inline int excluded_thermal_cell(const State* source,int i,int j,MeshView m,Parameters p) {
    if(p.shock_flags)return p.shock_flags[m.cell(i,j)]|p.shock_flags[m.cell(i-1,j)]|p.shock_flags[m.cell(i+1,j)]|p.shock_flags[m.cell(i,j-1)]|p.shock_flags[m.cell(i,j+1)];
    return entropy_exclusion_sensor(source,i,j,m,p)|entropy_exclusion_sensor(source,i-1,j,m,p)|entropy_exclusion_sensor(source,i+1,j,m,p)|entropy_exclusion_sensor(source,i,j-1,m,p)|entropy_exclusion_sensor(source,i,j+1,m,p);
}
enum ThermalFailureReason {
    THERMAL_OK=0,THERMAL_NONFINITE_STATE=1,THERMAL_NONPOSITIVE_DENSITY=2,
    THERMAL_NONPOSITIVE_AUXILIARY=3,THERMAL_PRESSURE_CONSTRUCTION=4,
    THERMAL_NONPOSITIVE_ENERGY_PRESSURE=5,THERMAL_FINAL_POSTCONDITION=6
};
HD inline bool thermal_recovery_failure(int* reason,int value) {
    if(reason)*reason=value;
    return false;
}
HD inline bool recover_thermal(int k,State* dest,const State* source,MeshView m,Parameters p,
                               double* floors,double* floor_energy,ThermalAccounting* accounting,double weight,int* failure_reason=nullptr) {
    int i=k%m.nx+NG,j=k/m.nx+NG,id=m.cell(i,j);State q=dest[id];
    if(failure_reason)*failure_reason=THERMAL_OK;
    for(int c=0;c<NV;++c)if(!finite_value(q[c]))return thermal_recovery_failure(failure_reason,THERMAL_NONFINITE_STATE);
    if(!(q[RHO]>0))return thermal_recovery_failure(failure_reason,THERMAL_NONPOSITIVE_DENSITY);
    double radius=m.radius(i,j,p),pressure=raw_pressure(q,radius,p);
    double internal=pressure/(p.gamma-1);
    bool use_entropy=false;
    if(p.entropy_switch==3&&p.internal_energy_sync==1) {
        // Enzo-E "modern": do not repeatedly seed a small, unreliable
        // E-K-M residual into the independently evolved gas energy merely
        // because a neighbouring pressure gradient was flagged as a shock.
        // All three conditions are required (including the half-aux bound).
        // This is a named dual-energy discretization, not a positivity proof
        // or an entropy equation guaranteed to capture every low-beta shock.
        use_entropy=!modern_energy_reliable(q,radius,p);
    } else use_entropy=p.entropy_switch&&!excluded_thermal_cell(source,i,j,m,p)
        &&(p.entropy_switch==1||internal<p.entropy_ratio*fabs(q[ENE]));
    if(use_entropy) {
        if(!(p.entropy_switch==3?q[EINT]>0:q[ENTR]>0))return thermal_recovery_failure(failure_reason,THERMAL_NONPOSITIVE_AUXILIARY);
        double ps=p.entropy_switch==3?(p.gamma-1)*q[EINT]:q[ENTR]*pow(q[RHO],p.gamma-1),old=q[ENE];
        if(!finite_value(ps)||!set_pressure(q,radius,fmax(ps,p.p_floor),p))return thermal_recovery_failure(failure_reason,THERMAL_PRESSURE_CONSTRUCTION);
        double floor_add=ps<p.p_floor?(p.p_floor-ps)/(p.gamma-1):0.;
        double correction=weight*(q[ENE]-old-floor_add)*m.volume(i,j,p);
        if(p.entropy_switch==3){accounting[k].internal_events+=1;accounting[k].internal_energy+=correction;}
        else {accounting[k].entropy_events+=1;accounting[k].entropy_energy+=correction;}
        pressure=raw_pressure(q,radius,p);dest[id]=q;
        if(ps<p.p_floor) {floors[k]+=1;floor_energy[k]+=weight*(p.p_floor-ps)/(p.gamma-1)*m.volume(i,j,p);q[ENTR]=pressure*pow(q[RHO],1-p.gamma);}
    } else {
        // Some reference solvers allow an explicitly nonconservative pressure
        // floor in energy-authoritative cells. Keep it selectable and fully
        // accounted; it must never hide nonfinite data or an invalid selected
        // auxiliary channel. The default remains rollback/retry for p<=0.
        if(!finite_value(pressure)||(!(pressure>0)&&!p.negative_pressure_floor))
            return thermal_recovery_failure(failure_reason,THERMAL_NONPOSITIVE_ENERGY_PRESSURE);
        if(pressure<p.p_floor) {
            post_ct_pressure_floor(k,dest,m,p,floors,floor_energy,weight);q=dest[id];pressure=raw_pressure(q,radius,p);
            if(!(pressure>0)||!finite_value(pressure))return thermal_recovery_failure(failure_reason,THERMAL_PRESSURE_CONSTRUCTION);
        }
    }
    // Legacy energy is authoritative in flagged shocks. Modern synchronization
    // instead uses its local thermal-reliability conditions.
    // Keep advected entropy in smooth entropy cells to avoid repeated reseeding.
    if(!use_entropy||p.entropy_switch==3)q[ENTR]=pressure*pow(q[RHO],1-p.gamma);
    if(!use_entropy||p.entropy_switch!=3)q[EINT]=pressure/(p.gamma-1);
    if(q[EINT]<p.p_floor/(p.gamma-1))q[EINT]=pressure/(p.gamma-1);
    dest[id]=q;
    for(int c=0;c<NV;++c)if(!finite_value(q[c]))return thermal_recovery_failure(failure_reason,THERMAL_NONFINITE_STATE);
    if(!(finite_value(raw_pressure(q,radius,p))&&raw_pressure(q,radius,p)>0&&q[ENTR]>0&&q[EINT]>0))
        return thermal_recovery_failure(failure_reason,THERMAL_FINAL_POSTCONDITION);
    return true;
}
struct StageFailure {
    int active=0,stage=0,i=0,j=0;
    double dt=0,pressure=0;
    State source,rhs,failed;
    int reason=THERMAL_OK,thermal_exclusion=0;
    double pre_ct_pressure=0,ct_delta_energy=0,auxiliary_pressure=0;
};
inline const char* thermal_failure_name(int reason) {
    switch(reason) {
        case THERMAL_NONFINITE_STATE:return "nonfinite_state";
        case THERMAL_NONPOSITIVE_DENSITY:return "nonpositive_density";
        case THERMAL_NONPOSITIVE_AUXILIARY:return "nonpositive_auxiliary";
        case THERMAL_PRESSURE_CONSTRUCTION:return "pressure_construction";
        case THERMAL_NONPOSITIVE_ENERGY_PRESSURE:return "nonpositive_energy_pressure";
        case THERMAL_FINAL_POSTCONDITION:return "final_postcondition";
        default:return "none";
    }
}
HD inline void capture_failure(StageFailure* record,int stage,int i,int j,double dt,const State& source,const State& rhs,const State& failed,double pressure,
                               int reason=THERMAL_OK,int thermal_exclusion=0,double pre_ct_pressure=0,double ct_delta_energy=0,double auxiliary_pressure=0) {
#ifdef __CUDA_ARCH__
    bool first=atomicCAS(&record->active,0,1)==0;
#else
    bool first=record->active==0;if(first)record->active=1;
#endif
    if(first){record->stage=stage;record->i=i;record->j=j;record->dt=dt;record->pressure=pressure;record->source=source;record->rhs=rhs;record->failed=failed;
        record->reason=reason;record->thermal_exclusion=thermal_exclusion;record->pre_ct_pressure=pre_ct_pressure;
        record->ct_delta_energy=ct_delta_energy;record->auxiliary_pressure=auxiliary_pressure;}
}
HD inline void rk_potential(int k,const double* base,const double* in,const double* rhs,double* out,MeshView m,double dt,double w0,double w1) {
    int i=k%(m.nx+1)+NG,j=k/(m.nx+1)+NG,id=m.vertex(i,j);out[id]=w0*base[id]+w1*(in[id]+dt*rhs[id]);
}
