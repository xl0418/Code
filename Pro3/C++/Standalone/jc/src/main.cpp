#include <iostream>
#include <stdexcept>
#include "jc.h"
#include "cmd_line.h"


const char* JCVersion = "jc version 1.0.2";
const char* JCReportBugs = "report bugs to h.hildenbrandt@rug.nl, please";


const char* JCHelp = R"(Usage: jc [OPTION]... [OPTIONAL PARAMETER]... PARAMETER... file
Options:
  --help           prints this text and exits
  --version        prints version information and exits
  -v, --verbose    event output to console
  -i, --implicit   run spacial implicit model
  -d, --dominant   run model at least until first specie isn't dominant

Optional parameter:
  L                grid size, defaults to 333
  v                speciation rate, defaults to 0.00001
  log              log interval [events] for M, D and R, defaults to 10
  log_first        always log the first n events, defaults to 0

Parameter:
  phi              strength of phylogenetic effect
  psi              strength of abundance effect
  sA               abundance Gaussian standard deviation
  sB               distribution Gaussian standard deviation
  ticks            time steps (turnovers)
  file             data file name

If option -i, --implicit is given, parameter sA and sB are ignored.

If option -d, --dominant is given, 'ticks' are the additional time steps
after first specie isn't dominant anymore.

Example:
jc phi=0.01 psi=0.1 sA=5 sB=0.5 ticks=10e7 seed=96 file="res.m"
)";


void JCModel(int argc, const char** argv, jc::Parameter& param)
{
  cmd::parse_optional_arg("L", param.L, argc, argv);
  cmd::parse_optional_arg("v", param.v, argc, argv);
  cmd::parse_optional_arg<int64_t, double>("log", param.log_interval, argc, argv);
  cmd::parse_optional_arg("log_first", param.log_first, argc, argv);
  cmd::parse_required_arg("phi", param.phi, argc, argv);
  cmd::parse_required_arg("psi", param.psi, argc, argv);
  if (param.implicit)
  {
    param.sigmaA = 100000.0;
    param.sigmaB = 100000.0;
  } 
  else 
  {
    cmd::parse_required_arg("sA", param.sigmaA, argc, argv);
    cmd::parse_required_arg("sB", param.sigmaB, argc, argv);
  }
  cmd::parse_required_arg<int64_t, double>("ticks", param.ticks, argc, argv);
  cmd::parse_required_arg("file", param.filename, argc, argv);
  if (param.L % 2 == 0) throw std::invalid_argument("L shall be odd");
  jc::RunModel(param);
}


void JCModel(int argc, const char** argv, jc::Parameter& param, bool /* continuation */)
{
  cmd::parse_required_arg<int64_t, double>("ticks", param.ticks, argc, argv);
  cmd::parse_required_arg("file", param.filename, argc, argv);
  throw std::invalid_argument("simulation continuation not yet implemented");
}


int main(int argc, const char** argv)
{
  // Default parameter, see declaration of struct Param
  jc::Parameter param;

  try
  {
    bool flag = false;
    cmd::parse_cmd_flag("--help", flag, argc, argv);
    if (flag)
    {
      std::cout << JCHelp << std::endl;
      return 0;
    }
    cmd::parse_cmd_flag("--version", flag, argc, argv);
    if (flag)
    {
      std::cout << JCVersion << '\n';
      std::cout << JCReportBugs << std::endl;
      return 0;
    }
    cmd::parse_cmd_flag("-v", param.verbose, argc, argv);
    cmd::parse_cmd_flag("--verbose", param.verbose, argc, argv);
    cmd::parse_cmd_flag("-c", param.continuation, argc, argv);
    cmd::parse_cmd_flag("--cont", param.continuation, argc, argv);
    cmd::parse_cmd_flag("-i", param.implicit, argc, argv);
    cmd::parse_cmd_flag("--implicit", param.implicit, argc, argv);
    cmd::parse_cmd_flag("-d", param.dominant, argc, argv);
    cmd::parse_cmd_flag("--dominant", param.dominant, argc, argv);
    cmd::parse_cmd_flag("--profile", param.profile, argc, argv);
    param.continuation ? JCModel(argc, argv, param, true) : JCModel(argc, argv, param);
    if (param.verbose) std::cout << "Regards\n";
    return 0;
  }
  catch (cmd::parse_error& e)
  {
    std::cerr << "jc: fatal error: " << e.what() << std::endl;
  }
  catch (std::exception& e)
  {
    std::cerr << "jc: fatal error: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "jc: fatal error: unknown exception" << std::endl;
  }
  return -1;
}
