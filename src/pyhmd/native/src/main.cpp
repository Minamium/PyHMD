#include "pyhmd/runtime/application.hpp"
#include "pyhmd/runtime/protocol.hpp"

int main(int argc,char** argv) {
    if(argc==2&&std::string(argv[1])=="--describe"){pyhmd::describe();return 0;}
    if(argc>=2&&std::string(argv[1])=="--validate-config")return pyhmd::validate_configuration(argc-1,argv+1);
    return pyhmd::run_application(argc,argv);
}
