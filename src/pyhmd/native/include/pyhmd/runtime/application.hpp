#pragma once
#include "pyhmd/io/output.hpp"

namespace pyhmd {
inline int run_application(int argc,char** argv) {
    std::unique_ptr<World> world;
    std::string active_output;double accepted_time=0;int accepted_step=0;
    try {
        world=std::make_unique<World>(argc,argv);Setup setup=parse_setup(argc,argv);Mesh mesh(setup.p,*world);HostMesh host(mesh);Solver solver(setup.p,mesh.view,*world);
        double time=0,dt=0,next_output=setup.p.output_interval;int step=0;
        if(setup.restart.empty())initialize(solver,host);
        else {
            const auto clock=pymhd_checkpoint::read(solver,setup.restart,setup.p.output_interval);
            time=clock.time;dt=clock.dt;next_output=clock.next_output;step=static_cast<int>(clock.step);
            if(time>=setup.p.t_end)throw std::runtime_error("restart time must be less than t_end");
            if(step>=setup.p.max_steps)throw std::runtime_error("restart step must be less than max_steps");
        }
        Output output(setup,*world,mesh,host);solver.failure_directory=setup.output;
        active_output=setup.output;accepted_time=time;accepted_step=step;
        const auto start=std::chrono::steady_clock::now();auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
        double last_accepted_dt=dt;
        auto clock=[&](){return pymhd_checkpoint::Clock{time,last_accepted_dt,next_output,setup.p.output_interval,step};};
        bool saved_first=false,saved_severe=false,saved_limit=false;
        if(setup.checkpoint_on_retry)solver.rejected_trial_hook=[&](const Solver& rejected,const StageFailure& failure,int attempt){
            std::vector<const char*> names;
            if(!saved_first){saved_first=true;names.push_back("retry_first");}
            else if(attempt>=7&&!saved_limit){saved_limit=true;names.push_back("retry_limit");}
            else if(attempt>=3&&!saved_severe){saved_severe=true;names.push_back("retry_severe");}
            if(attempt==rejected.p.max_step_retries)names.push_back("retry_exhausted");
            for(const char* name:names){pymhd_checkpoint::write(rejected,fs::path(setup.output)/name,clock(),true,&failure,attempt);if(world->rank==0)std::cerr<<"Saved last accepted state before rejection: "<<(fs::path(setup.output)/name)<<'\n';}
        };
        double next_checkpoint=setup.checkpoint_interval>0?(std::floor(time/setup.checkpoint_interval)+1)*setup.checkpoint_interval:INFINITY;
        int retry_streak=0;std::string stop_reason="max_steps";
        output.write(solver,step,time,dt,elapsed(),true);
        while(time<setup.p.t_end && step<setup.p.max_steps) {
            const double cfl_dt=solver.timestep();dt=std::min(cfl_dt,setup.p.t_end-time);
            if(next_output>time)dt=std::min(dt,next_output-time);
            if(!(dt>1e-16)||!std::isfinite(dt))throw std::runtime_error("time step collapsed");
            const double requested_dt=dt;dt=solver.step(dt);last_accepted_dt=dt;time+=dt;++step;
            double clock_tolerance=8*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(time));
            if(std::abs(time-next_output)<=clock_tolerance)time=next_output;
            if(std::abs(time-setup.p.t_end)<=clock_tolerance)time=setup.p.t_end;
            accepted_time=time;accepted_step=step;
            bool snapshot=time>=next_output-1e-12*std::max(1.,time)||time>=setup.p.t_end;
            retry_streak=solver.last_step_retries>=setup.retry_streak_threshold?retry_streak+1:0;
            const bool retry_stop=setup.retry_streak_limit>0&&retry_streak>=setup.retry_streak_limit;
            if(snapshot||step%setup.p.log_interval==0||step==setup.p.max_steps||retry_stop) {
                output.write(solver,step,time,dt,elapsed(),snapshot||step==setup.p.max_steps||retry_stop,cfl_dt,requested_dt);
                // A final snapshot between scheduled outputs must not skip the
                // next output time when this state is subsequently restarted.
                while(time>=next_output-1e-12*std::max(1.,time))next_output+=setup.p.output_interval;
            }
            if(time>=next_checkpoint||retry_stop||time>=setup.p.t_end||step>=setup.p.max_steps) {
                const auto target=retry_stop?fs::path(setup.output)/"retry_stop":fs::path(setup.output)/("checkpoint_step"+std::to_string(step));
                if(setup.checkpoint_interval>0||retry_stop)pymhd_checkpoint::write(solver,target,clock());
                if(setup.checkpoint_interval>0)next_checkpoint=(std::floor(time/setup.checkpoint_interval)+1)*setup.checkpoint_interval;
            }
            if(retry_stop){stop_reason="persistent_retries";if(world->rank==0)std::cerr<<"Stopping after "<<retry_streak<<" consecutive steps requiring at least "<<setup.retry_streak_threshold<<" rejected trials; see retry_stop checkpoint.\n";break;}
        }
        if(time>=setup.p.t_end)stop_reason="t_end";
        if(world->rank==0) {
            std::ofstream status(fs::path(setup.output)/"status.json");status<<std::setprecision(17)<<"{\"complete\":"<<(time>=setup.p.t_end?"true":"false")<<",\"reason\":\""<<stop_reason<<"\",\"time\":"<<time<<",\"step\":"<<step<<",\"wall_seconds\":"<<elapsed()<<"}\n";
        }
        return time>=setup.p.t_end?0:stop_reason=="persistent_retries"?4:3;
    } catch(const std::exception& e) {
        if(world&&world->rank==0&&!active_output.empty()) {
            // The exception path must not leave an incomplete run without a
            // machine-readable status. Do not overwrite another run when
            // parsing/setup or output-directory ownership checks failed.
            std::ofstream status(fs::path(active_output)/"status.json");
            status<<std::setprecision(17)<<"{\"complete\":false,\"reason\":\"exception\",\"time\":"<<accepted_time<<",\"step\":"<<accepted_step<<"}\n";
            status.close();
        }
        std::cerr<<"PyMHD error: "<<e.what()<<'\n';if(world)world->abort(2);return 2;
    }
}

} // namespace pyhmd
