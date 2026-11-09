#pragma once
#include "pyhmd/problems/builtin.hpp"
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

struct Setup {
    Parameters p;
    std::string output="runs/thin_disk",config_path,restart;
    double checkpoint_interval=0;
    int checkpoint_on_retry=1;
    int midplane_diagnostics=0;
    int retry_streak_limit=64,retry_streak_threshold=4;
};

inline std::string setup_trim(const std::string& text) {
    const auto first=text.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return {};
    return text.substr(first,text.find_last_not_of(" \t\r\n")-first+1);
}
inline std::string setup_lower(std::string value) {
    for(char& c:value)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return value;
}
inline double setup_double(const std::string& value,const std::string& key) {
    size_t used=0;double number;
    try {number=std::stod(value,&used);}catch(const std::exception&){throw std::runtime_error("invalid number for "+key+": "+value);}
    if(used!=value.size() || !std::isfinite(number))throw std::runtime_error("invalid finite number for "+key+": "+value);
    return number;
}
inline int setup_integer(const std::string& value,const std::string& key) {
    size_t used=0;long long number;
    try {number=std::stoll(value,&used);}catch(const std::exception&){throw std::runtime_error("invalid integer for "+key+": "+value);}
    if(used!=value.size() || number<std::numeric_limits<int>::min() || number>std::numeric_limits<int>::max())throw std::runtime_error("invalid integer for "+key+": "+value);
    return static_cast<int>(number);
}
inline int setup_boolean(const std::string& value,const std::string& key) {
    const auto name=setup_lower(value);
    if(name=="1" || name=="true" || name=="yes" || name=="on")return 1;
    if(name=="0" || name=="false" || name=="no" || name=="off" || name=="none")return 0;
    throw std::runtime_error("invalid boolean for "+key+": "+value);
}
inline int setup_case(const std::string& name) {
    const auto value=setup_lower(name);
    if(value=="uniform")return UNIFORM;
    if(value=="alfven")return ALFVEN;
    if(value=="diffusion")return DIFFUSION;
    if(value=="brio_wu" || value=="brio-wu")return BRIO_WU;
    if(value=="thin_disk" || value=="disk")return DISK;
    if(value=="oblique_alfven")return OBLIQUE_ALFVEN;
    if(value=="linear_wave")return LINEAR_WAVE;
    if(value=="field_loop")return FIELD_LOOP;
    if(value=="rj2a")return RJ2A;
    if(value=="orszag_tang")return ORSZAG_TANG;
    if(value=="blast_low_beta")return BLAST_LOW_BETA;
    throw std::runtime_error("unknown case: "+name);
}
inline int setup_boundary(const std::string& name) {
    const auto value=setup_lower(name);
    if(value=="outflow")return OUTFLOW;
    if(value=="periodic")return PERIODIC;
    if(value=="reflect" || value=="reflecting")return REFLECT;
    if(value=="axis")return AXIS;
    if(value=="equator")return EQUATOR;
    if(value=="paper_inner")return PAPER_INNER;
    if(value=="paper_outer")return PAPER_OUTER;
    throw std::runtime_error("unknown boundary: "+name);
}

inline Parameters setup_case_defaults(int problem) {
    Parameters p;p.problem=problem;
    if(problem==DISK) {
        p.nx=600;p.ny=128;p.spherical=1;p.log_r=1;
        p.xlo=1;p.xhi=1500;p.ylo=0;p.yhi=PI/2;
        p.lower1=PAPER_INNER;p.upper1=PAPER_OUTER;p.lower2=AXIS;p.upper2=EQUATOR;
        p.gamma=5./3;p.cfl=.4;p.gm=1;p.t_end=100;p.output_interval=10;
        p.eta_model=1;p.cooling=1;p.h=.1;p.corona=.001;
        p.entropy_switch=3;p.ct_energy_correction=1;
        p.internal_energy_sync=1;
        p.ct_emf_scheme=1;
        p.b0=.1*sqrt(.02);p.alpha_m=1.65;p.chi=.5;p.mu0=.01;p.mu_ref=.01;p.bend=.5;
        return p;
    }
    p.gm=0;p.eta=0;p.eta_model=0;p.cooling=0;p.log_r=0;
    p.lower1=p.upper1=p.lower2=p.upper2=OUTFLOW;
    p.t_end=.2;p.output_interval=.05;p.log_interval=50;
    if(problem==UNIFORM) {
        p.nx=64;p.ny=32;p.spherical=1;
        p.xlo=1;p.xhi=2;p.ylo=.2;p.yhi=PI-.2;
    }else {
        p.spherical=0;p.nx=128;p.ny=4;p.xlo=0;p.xhi=1;p.ylo=0;p.yhi=.125;
        p.lower2=p.upper2=PERIODIC;
        if(problem==ALFVEN || problem==DIFFUSION)p.lower1=p.upper1=PERIODIC;
        if(problem==ALFVEN)p.wave_amplitude=.1;
        if(problem==DIFFUSION){p.wave_amplitude=.001;p.eta=.01;}
        if(problem==BRIO_WU){p.gamma=2;p.nx=256;p.t_end=.1;p.output_interval=.025;}
    }
    if(benchmark_problem(problem)) {
        p.spherical=p.log_r=0;p.nx=64;p.ny=32;p.gamma=5./3;p.t_end=.1;
        p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;p.ct_emf_scheme=1;
        p.xlo=p.ylo=0;p.xhi=sqrt(5.);p.yhi=sqrt(5.)/2;
        if(problem==OBLIQUE_ALFVEN){p.wave_amplitude=.1;p.t_end=1;}
        if(problem==LINEAR_WAVE){p.wave_amplitude=1e-6;p.t_end=.5;}
        if(problem==FIELD_LOOP){p.xlo=-1;p.xhi=1;p.ylo=-.5;p.yhi=.5;p.wave_amplitude=.001;p.t_end=1;}
        if(problem==ORSZAG_TANG||problem==BLAST_LOW_BETA){p.xlo=p.ylo=-.5;p.xhi=p.yhi=.5;p.ny=p.nx;p.t_end=.1;}
        if(problem==RJ2A){p.xlo=-.5;p.xhi=.5;p.ylo=0;p.yhi=.125;p.nx=128;p.ny=4;p.lower1=p.upper1=OUTFLOW;p.t_end=.2;}
    }
    return p;
}

inline std::map<std::string,std::string> setup_read_ini(const std::string& path) {
    std::ifstream input(path);
    if(!input)throw std::runtime_error("cannot open configuration: "+path);
    std::map<std::string,std::string> settings;std::string line;int line_number=0;
    while(std::getline(input,line)) {
        ++line_number;
        const auto comment=line.find_first_of("#;");
        if(comment!=std::string::npos)line.resize(comment);
        line=setup_trim(line);
        if(line.empty())continue;
        if(line.front()=='[') {
            if(line.back()!=']')throw std::runtime_error(path+":"+std::to_string(line_number)+": malformed section");
            continue;
        }
        const auto equal=line.find('=');
        if(equal==std::string::npos)throw std::runtime_error(path+":"+std::to_string(line_number)+": expected key=value");
        const std::string key=setup_lower(setup_trim(line.substr(0,equal)));
        std::string value=setup_trim(line.substr(equal+1));
        if(key.empty() || value.empty())throw std::runtime_error(path+":"+std::to_string(line_number)+": empty key/value");
        if(value.size()>=2 && ((value.front()=='"' && value.back()=='"') || (value.front()=='\'' && value.back()=='\'')))value=value.substr(1,value.size()-2);
        if(!settings.emplace(key,value).second)throw std::runtime_error(path+":"+std::to_string(line_number)+": duplicate key "+key);
    }
    return settings;
}

inline void setup_assign(Setup& setup,const std::string& key,const std::string& value) {
    Parameters& p=setup.p;
    if(key=="case"){p.problem=setup_case(value);return;}
    if(key=="wave_mode") {
        const auto mode=setup_lower(value);
        if(mode=="fast")p.wave_mode=0;else if(mode=="alfven")p.wave_mode=1;
        else if(mode=="slow")p.wave_mode=2;else if(mode=="entropy")p.wave_mode=3;
        else throw std::runtime_error("unknown linear wave mode: "+value);
        return;
    }
    if(key=="noise_seed"){p.noise_seed=setup_integer(value,key);return;}
    if(key=="noise_amplitude"){p.noise_amplitude=setup_double(value,key);return;}
    if(key=="output"){setup.output=value;return;}
    if(key=="restart"){setup.restart=value;return;}
    if(key=="checkpoint_interval"){setup.checkpoint_interval=setup_double(value,key);return;}
    if(key=="checkpoint_on_retry"){setup.checkpoint_on_retry=setup_boolean(value,key);return;}
    if(key=="midplane_diagnostics"){setup.midplane_diagnostics=setup_boolean(value,key);return;}
    if(key=="transport_diagnostics"){p.transport_diagnostics=setup_boolean(value,key);return;}
    if(key=="retry_streak_limit"){setup.retry_streak_limit=setup_integer(value,key);return;}
    if(key=="retry_streak_threshold"){setup.retry_streak_threshold=setup_integer(value,key);return;}
    if(key=="internal_energy_sync") {
        auto name=setup_lower(value);
        if(name=="modern")p.internal_energy_sync=1;
        else if(name=="legacy")p.internal_energy_sync=0;
        else throw std::runtime_error("internal_energy_sync must be modern or legacy");
        return;
    }
    if(key=="dual_energy_eta"){p.dual_energy_eta=setup_double(value,key);return;}
    if(key=="cfl_staggered"){p.cfl_staggered=setup_boolean(value,key);return;}
    if(key=="inner_normal_field") {
        auto name=setup_lower(value);
        if(name=="evolve")p.inner_flux_project=0;
        else if(name=="project")p.inner_flux_project=1;
        else throw std::runtime_error("inner_normal_field must be evolve or project");
        return;
    }
    if(key=="ct_emf_scheme") {
        auto name=setup_lower(value);
        if(name=="corner_hll")p.ct_emf_scheme=0;
        else if(name=="face_hll")p.ct_emf_scheme=1;
        else throw std::runtime_error("ct_emf_scheme must be corner_hll or face_hll");
        return;
    }
    if(key=="geometry") {
        const auto name=setup_lower(value);
        if(name=="spherical")p.spherical=1;
        else if(name=="cartesian")p.spherical=0;
        else throw std::runtime_error("unknown geometry: "+value);
        return;
    }
    if(key=="radial_grid") {
        const auto name=setup_lower(value);
        if(name=="log" || name=="logarithmic")p.log_r=1;
        else if(name=="uniform")p.log_r=0;
        else throw std::runtime_error("unknown radial_grid: "+value);
        return;
    }
    if(key=="log_r"){p.log_r=setup_boolean(value,key);return;}
    if(key=="cooling"){p.cooling=setup_boolean(value,key);return;}
    if(key=="ct_energy_correction"){p.ct_energy_correction=setup_boolean(value,key);return;}
    if(key=="reconstruction_cache"){p.reconstruction_cache=setup_boolean(value,key);return;}
    if(key=="metric_ppm"){p.metric_ppm=setup_boolean(value,key);return;}
    if(key=="retry_flat_same_dt"){p.retry_flat_same_dt=setup_boolean(value,key);return;}
    if(key=="negative_pressure_policy") {
        const auto name=setup_lower(value);
        if(name=="retry"||name=="strict"||name=="0")p.negative_pressure_floor=0;
        else if(name=="floor"||name=="1")p.negative_pressure_floor=1;
        else throw std::runtime_error("unknown negative_pressure_policy: "+value);
        return;
    }
    if(key=="initial_equilibrium") {
        const auto name=setup_lower(value);
        if(name=="hydrostatic"||name=="0")p.initial_equilibrium=0;
        else if(name=="magnetohydrostatic"||name=="lorentz_balanced"||name=="1")p.initial_equilibrium=1;
        else throw std::runtime_error("unknown initial_equilibrium: "+value);
        return;
    }
    if(key=="entropy_switch"||key=="thermal_recovery") {
        const auto name=setup_lower(value);
        if(name=="off"||name=="no"||name=="0")p.entropy_switch=0;
        else if(name=="selective"||name=="yes"||name=="1")p.entropy_switch=1;
        else if(name=="dual_energy"||name=="2")p.entropy_switch=2;
        else if(name=="internal_energy"||name=="3")p.entropy_switch=3;
        else throw std::runtime_error("unknown entropy_switch: "+value);
        return;
    }
    if(key=="eta_model") {
        const auto name=setup_lower(value);
        if(name=="constant" || name=="0")p.eta_model=0;
        else if(name=="sf2014" || name=="paper" || name=="1")p.eta_model=1;
        else if(name=="strong" || name=="2")p.eta_model=2;
        else throw std::runtime_error("unknown eta_model: "+value);
        return;
    }
    if(key=="lower1"){p.lower1=setup_boundary(value);return;}
    if(key=="upper1"){p.upper1=setup_boundary(value);return;}
    if(key=="lower2"){p.lower2=setup_boundary(value);return;}
    if(key=="upper2"){p.upper2=setup_boundary(value);return;}
    const std::map<std::string,int Parameters::*> integers={
        {"nx1",&Parameters::nx},{"nx2",&Parameters::ny},{"max_steps",&Parameters::max_steps},{"log_interval",&Parameters::log_interval},
        {"max_step_retries",&Parameters::max_step_retries}};
    auto integer=integers.find(key);
    if(integer!=integers.end()){p.*(integer->second)=setup_integer(value,key);return;}
    const std::map<std::string,double Parameters::*> numbers={
        {"x1_min",&Parameters::xlo},{"x1_max",&Parameters::xhi},{"x2_min",&Parameters::ylo},{"x2_max",&Parameters::yhi},
        {"gamma",&Parameters::gamma},{"cfl",&Parameters::cfl},{"t_end",&Parameters::t_end},{"output_interval",&Parameters::output_interval},
        {"gm",&Parameters::gm},{"eta",&Parameters::eta},{"rho_floor",&Parameters::rho_floor},{"p_floor",&Parameters::p_floor},
        {"h",&Parameters::h},{"corona",&Parameters::corona},{"b0",&Parameters::b0},{"mu0",&Parameters::mu0},
        {"bend",&Parameters::bend},{"alpha_m",&Parameters::alpha_m},{"chi",&Parameters::chi},{"mu_ref",&Parameters::mu_ref},
        {"wave_amplitude",&Parameters::wave_amplitude},{"entropy_ratio",&Parameters::entropy_ratio},
        {"entropy_shock_threshold",&Parameters::entropy_shock_threshold},{"entropy_contact_threshold",&Parameters::entropy_contact_threshold}};
    auto number=numbers.find(key);
    if(number!=numbers.end()){p.*(number->second)=setup_double(value,key);return;}
    throw std::runtime_error("unknown configuration key: "+key);
}

inline void setup_validate(const Setup& setup) {
    const Parameters& p=setup.p;
    auto require=[](bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);};
    require(p.nx>=NG && p.ny>=NG,"nx1 and nx2 must each be at least the halo width (3)");
    require(p.xhi>p.xlo && p.yhi>p.ylo,"grid upper bounds must exceed lower bounds");
    require(p.gamma>1,"gamma must exceed 1");
    require(p.cfl>0 && p.cfl<=.5,"cfl must satisfy 0 < cfl <= 0.5 for this unsplit 2D RK scheme");
    require(p.t_end>0 && p.max_steps>0,"t_end and max_steps must be positive");
    require(p.output_interval>0 && p.log_interval>0,"output_interval and log_interval must be positive");
    require(p.rho_floor>0 && p.p_floor>0,"density and pressure floors must be positive");
    require(p.entropy_ratio>0 && p.entropy_ratio<1 && p.entropy_shock_threshold>0,"invalid entropy switch thresholds");
    require(p.dual_energy_eta>0&&p.dual_energy_eta<1,"dual_energy_eta must lie between zero and one");
    require(p.entropy_contact_threshold>0 && p.entropy_contact_threshold<1,"entropy_contact_threshold must lie between zero and one");
    require(p.max_step_retries>=0 && p.max_step_retries<=20,"max_step_retries must be between 0 and 20");
    require(p.gm>=0 && p.eta>=0,"gm and eta must be nonnegative");
    require(p.h>0 && p.h<1,"disk h must satisfy 0 < h < 1");
    require(p.corona>0 && p.alpha_m>=0 && p.chi>=0 && p.mu_ref>0 && p.mu0>=0 && p.bend>0,"invalid disk/corona/resistivity parameter");
    require(p.wave_amplitude>=0,"wave_amplitude must be nonnegative");
    require(p.noise_seed>=0&&p.noise_amplitude>=0&&p.noise_amplitude<=1e-6,"invalid initial noise controls");
    require(p.noise_amplitude==0||p.problem==DISK,"initial noise is supported only for thin_disk");
    require(p.wave_mode==0||p.problem==LINEAR_WAVE,"wave_mode applies only to linear_wave");
    require(p.problem!=LINEAR_WAVE||(p.wave_amplitude>0&&p.wave_amplitude<=1e-3),"linear_wave requires 0 < wave_amplitude <= 1e-3");
    if(benchmark_problem(p.problem)) {
        require(!p.spherical&&!p.log_r&&p.eta_model==0&&p.eta==0&&p.gm==0&&!p.cooling,"benchmarks require Cartesian uniform ideal MHD without gravity/cooling");
        require(std::abs(p.gamma-5./3)<1e-14,"benchmarks use gamma=5/3");
        if(p.problem!=RJ2A)require(p.lower1==PERIODIC&&p.upper1==PERIODIC&&p.lower2==PERIODIC&&p.upper2==PERIODIC,"multidimensional benchmarks require periodic boundaries");
    }
    require(!setup.output.empty(),"output directory must not be empty");
    require(std::isfinite(setup.checkpoint_interval)&&setup.checkpoint_interval>=0,"checkpoint_interval must be finite and nonnegative");
    require(!p.transport_diagnostics||p.ct_emf_scheme==1,"transport_diagnostics requires face_hll CT");
    require(!p.transport_diagnostics||!p.inner_flux_project,"transport_diagnostics requires inner_normal_field=evolve");
    require(!p.transport_diagnostics||(p.problem==DISK&&p.spherical&&p.log_r&&p.lower2==AXIS&&p.upper2==EQUATOR),
            "transport_diagnostics output requires the spherical logarithmic disk with axis/equator boundaries");
    require(setup.retry_streak_limit>=0&&setup.retry_streak_threshold>=1&&setup.retry_streak_threshold<=20,"invalid retry streak stop limits");
    require(!p.log_r || p.xlo>0,"logarithmic grid requires x1_min > 0");
    require((p.lower1==PERIODIC)==(p.upper1==PERIODIC),"x1 periodic boundaries must be paired");
    require((p.lower2==PERIODIC)==(p.upper2==PERIODIC),"x2 periodic boundaries must be paired");
    if(p.spherical) {
        require(p.xlo>0,"spherical grid requires a strictly positive inner radius");
        require(p.ylo>=0 && p.yhi<=PI,"spherical theta bounds must lie within [0,pi]");
        if(p.lower2==AXIS)require(std::abs(p.ylo)<1e-13,"lower axis boundary requires theta_min=0");
        if(p.upper2==AXIS)require(std::abs(p.yhi-PI)<1e-13,"upper axis boundary requires theta_max=pi");
        if(p.lower2==EQUATOR)require(std::abs(p.ylo-PI/2)<1e-13,"lower equator boundary requires theta_min=pi/2");
        if(p.upper2==EQUATOR)require(std::abs(p.yhi-PI/2)<1e-13,"upper equator boundary requires theta_max=pi/2");
    }else {
        require(p.lower2!=AXIS && p.upper2!=AXIS && p.lower2!=EQUATOR && p.upper2!=EQUATOR,"axis/equator boundaries require spherical geometry");
    }
    require(p.problem!=DISK || p.spherical,"thin_disk requires spherical geometry");
    require(!setup.midplane_diagnostics||(p.problem==DISK&&p.spherical&&p.upper2==EQUATOR&&p.xlo<1.5&&p.xhi>1.1),"midplane diagnostics require an equatorial disk intersecting r=[1.1,1.5]");
    require(p.initial_equilibrium==0||p.initial_equilibrium==1,"invalid initial equilibrium");
    require(p.eta_model==0 || (p.problem==DISK && p.spherical && p.upper2==EQUATOR && std::abs(p.yhi-PI/2)<1e-13),"SF2014 resistivity requires the upper-hemisphere thin_disk case with equatorial upper boundary");
    require(p.lower2!=EQUATOR && p.upper2!=AXIS,"lower-equator and upper-axis orientations are not implemented");
}

inline Setup parse_setup(int argc,char** argv) {
    std::map<std::string,std::string> command;
    std::string config_path;
    const std::map<std::string,std::string> options={{"--output","output"},{"--case","case"},{"--nx1","nx1"},{"--nx2","nx2"},
        {"--t-end","t_end"},{"--max-steps","max_steps"},{"--eta","eta"},{"--cfl","cfl"},
        {"--metric-ppm","metric_ppm"},{"--initial-equilibrium","initial_equilibrium"},
        {"--retry-flat-same-dt","retry_flat_same_dt"},
        {"--negative-pressure-policy","negative_pressure_policy"},
        {"--internal-energy-sync","internal_energy_sync"},
        {"--dual-energy-eta","dual_energy_eta"},
        {"--ct-emf-scheme","ct_emf_scheme"},
        {"--cfl-staggered","cfl_staggered"},
        {"--inner-normal-field","inner_normal_field"},
        {"--transport-diagnostics","transport_diagnostics"},
        {"--midplane-diagnostics","midplane_diagnostics"},
        {"--output-interval","output_interval"},{"--log-interval","log_interval"},
        {"--restart","restart"},{"--checkpoint-interval","checkpoint_interval"},{"--checkpoint-on-retry","checkpoint_on_retry"},
        {"--retry-streak-limit","retry_streak_limit"},{"--retry-streak-threshold","retry_streak_threshold"}};
    for(int index=1;index<argc;++index) {
        const std::string option=argv[index];
        if(option=="--help" || option=="-h") {
            std::cout<<"PyMHD: finite-volume PPM/HLL/SSP-RK3 resistive MHD with constrained transport\n"
                     <<"Usage: pymhd_[cpu|cuda] [--config FILE] [--case NAME] [--output DIR]\n"
                     <<"       [--nx1 N] [--nx2 N] [--t-end T] [--max-steps N] [--eta VALUE]\n"
                     <<"       [--cfl C] [--output-interval T] [--log-interval N]\n"
                     <<"       [--metric-ppm yes|no] [--initial-equilibrium hydrostatic|magnetohydrostatic]\n"
                     <<"       [--retry-flat-same-dt yes|no]\n"
                     <<"       [--negative-pressure-policy retry|floor]\n"
                     <<"       [--internal-energy-sync modern|legacy]\n"
                     <<"       [--dual-energy-eta VALUE]\n"
                     <<"       [--ct-emf-scheme corner_hll|face_hll]\n"
                     <<"       [--cfl-staggered yes|no]\n"
                     <<"       [--inner-normal-field evolve|project]\n"
                     <<"       [--transport-diagnostics yes|no]\n"
                     <<"       [--restart CHECKPOINT_DIR] [--checkpoint-interval T] [--checkpoint-on-retry yes|no]\n"
                     <<"       [--retry-streak-limit N] [--retry-streak-threshold N]\n"
                     <<"Cases: thin_disk, uniform, alfven, diffusion, brio_wu\n"
                     <<"       oblique_alfven, linear_wave, field_loop, rj2a, orszag_tang, blast_low_beta\n"
                     <<"Precedence: case defaults, then INI keys, then explicit CLI values.\n"
                     <<"INI sections are cosmetic; keys are globally unique and unknown keys fail.\n";
            std::exit(0);
        }
        if(option!="--config" && options.find(option)==options.end())throw std::runtime_error("unknown command-line option: "+option);
        if(index+1>=argc)throw std::runtime_error("missing value after "+option);
        const std::string value=argv[++index];
        if(option=="--config") {
            if(!config_path.empty())throw std::runtime_error("--config specified more than once");
            config_path=value;
        }else if(!command.emplace(options.at(option),value).second)throw std::runtime_error("duplicate command-line option: "+option);
    }
    std::map<std::string,std::string> settings;
    if(!config_path.empty())settings=setup_read_ini(config_path);
    int problem=DISK;
    if(settings.count("case"))problem=setup_case(settings.at("case"));
    if(command.count("case"))problem=setup_case(command.at("case"));
    Setup setup;setup.p=setup_case_defaults(problem);setup.config_path=config_path;
    for(const auto& item:settings)if(item.first!="case")setup_assign(setup,item.first,item.second);
    for(const auto& item:command)setup_assign(setup,item.first,item.second);
    if(setup.p.problem==DISK) {
        auto& p=setup.p;const bool explicit_b0=settings.count("b0");
        if(explicit_b0 && !settings.count("mu0"))p.mu0=p.b0*p.b0/(2*p.h*p.h*p.gm);
        const double consistent=p.h*std::sqrt(2*p.mu0*p.gm);
        if(explicit_b0 && std::abs(p.b0-consistent)>1e-12*std::max(1.,std::abs(consistent)))
            throw std::runtime_error("inconsistent b0 and mu0: require b0=h*sqrt(2*mu0*GM), or omit b0");
        p.b0=consistent;
    }
    setup_validate(setup);
    return setup;
}
