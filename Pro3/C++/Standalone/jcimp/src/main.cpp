#include <iostream>
#include <stdexcept>
#include "jc_base.h"
#include "cmd_line.h"


const char* JCVersion = "jc version 1.0.2";
const char* JCReportBugs = "report bugs to h.hildenbrandt@rug.nl, please";


const char* JCHelp = R"(Usage: jc [OPTION]... [OPTIONAL PARAMETER]... PARAMETER... file
Options:
  --help           prints this text and exits
  --version        prints version information and exits
  -v, --verbose    event output to console
  -i, --implicit   run spacial implicit model

Optional parameter:
  N                individuals, defaults to 100000
  v                speciation rate, defaults to 0.00001

Parameter:
  phi              strength of phylogenetic effect
  psi              strength of abundance effect
  sA               abundance Gaussian standard deviation
  sB               distribution Gaussian standard deviation
  ticks            time steps (turnovers)
  log              log interval [events_] for M, D and R
  file             data file name

If option -i, --implicit is given, parameter sA and sB are ignored.

Example:
jc phi=0.01 psi=0.1 sA=5 sB=0.5 ticks=1e7 file="res.m"
)";


void JCModel(int argc, const char** argv, jc::Parameter& param)
{
  cmd::parse_optional_arg<int, double>("N", param.N, argc, argv);
  cmd::parse_optional_arg("v", param.v, argc, argv);
  cmd::parse_required_arg<int64_t, double>("log", param.log_interval, argc, argv);
  cmd::parse_required_arg("phi", param.phi, argc, argv);
  cmd::parse_required_arg("psi", param.psi, argc, argv);
  cmd::parse_required_arg<int64_t, double>("ticks", param.ticks, argc, argv);
  cmd::parse_required_arg("file", param.filename, argc, argv);
  cmd::parse_required_arg("rep", param.rep, argc, argv);
  auto jcm = jc::CreateModel(param);
  jcm->run();
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
    cmd::parse_cmd_flag("-i", param.implicit, argc, argv);
    cmd::parse_cmd_flag("--implicit", param.implicit, argc, argv);
    cmd::parse_cmd_flag("-v", param.verbose, argc, argv);
    cmd::parse_cmd_flag("--verbose", param.verbose, argc, argv);
    JCModel(argc, argv, param);
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
