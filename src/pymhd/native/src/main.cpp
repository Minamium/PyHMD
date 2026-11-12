#include "pymhd/runtime/application.hpp"
#include "pymhd/runtime/protocol.hpp"

int main(int argc,char** argv) {
    if(argc==2&&std::string(argv[1])=="--describe"){pymhd::describe();return 0;}
    if(argc>=2&&std::string(argv[1])=="--validate-config")return pymhd::validate_configuration(argc-1,argv+1);
    return pymhd::run_application(argc,argv);
}
