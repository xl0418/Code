#include <stdexcept>
#include "rndutils.hpp"
#include <ctime>
#include <string>
#include "iomatlab.h"


extern const char* JCVersion;


namespace jc {


  // to facilitate parsing, these labels are inserted into the result stream.
  const char cpp_skip_label_read_back[] = "cpp_skip_label_read_back";
  const char cpp_skip_label_last_log[] = "cpp_skip_label_last_log";
  const char cpp_skip_label_event_log[] = "cpp_skip_label_event_log";
  const char cpp_skip_label_state[] = "cpp_skip_label_state";


  const char EventHeader[] = R"m(
%
% Speciation event log
%
% events(:,1)    The time the event occurred.
% events(:,2)    Number of species after the event.
% events(:,3:4)  Position on the grid.
% events(:,5)    New or extinct species.
% events(:,6)    Ancestor in case of speciation event, -1 otherwise.
%
% extinctions    Extinction events from events
% speciations    Speciation events from events
%
)m";


  std::ostream& operator<<(std::ostream& os, event const& e)
  {
    os << "  " << e.T << ' ' << e.NS << ' ' << e.pos.first << ' ' << e.pos.second << ' ' << e.sp << ' ' << e.ancestor << ';';
    return os;
  }


  std::istream& operator>>(std::istream& is, event& e)
  {
    is >> e.T >> e.NS >> e.pos.first >> e.pos.second >> e.sp >> e.ancestor;
    char delim; is >> delim;    // skip ';'
    return is;
  }


  std::ostream& operator<<(std::ostream& os, Parameter const& p)
  {
    os << p.L << ' ' << p.v << ' ';
    os << p.phi << ' ' << p.psi << ' ';
    os << p.sigmaA << ' ' << p.sigmaB << ' ';
    os << p.ticks << ' ';
    os << '"' << p.filename.c_str() << "\" ";
    os << p.log_interval << ' ' << p.elapsed_time;
    return os;
  }


  std::istream& operator<<(std::istream& is, Parameter& p)
  {
    is >> p.L >> p.v;
    is >> p.phi >> p.psi;
    is >> p.sigmaA >> p.sigmaB;
    is >> p.ticks;
    is >> p.filename;
    is >> p.log_interval >>p.elapsed_time;
    return is;
  }


  template <typename T>
  std::ostream& operator<<(std::ostream& os, square_buffer<T> const& s)
  {
    for (size_t i = 0; i < s.n(); ++i)
    {
      for (size_t j = 0; j < s.n(); ++j)
      {
        os << s(j, i) << ' ';
      }
      os << ";\n";
    }
    return os;
  }
  
  
  template <typename T>
  std::istream& operator>>(std::istream& is, square_buffer<T>& s)
  {
    for (size_t i = 0; i < s.n(); ++i)
    {
      for (size_t j = 0; j < s.n(); ++j)
      {
        is >> s(j, i);
      }
      char delim;
      is >> delim;    // skip ';'
    }
    return is;
  }


  MatLogger::MatLogger(Parameter const& param) : os_(param.filename)
  {
    if (!os_) throw std::invalid_argument((std::string("can't create data file ") + param.filename).c_str());

    auto today = std::time(nullptr);
    // spits out header
    os_ << "% " << param.filename.c_str();
    os_ << "\n% JanzenConnell result file.\n";
    os_ << "% " << JCVersion << '\n';
    os_ << "% Generated at " << std::ctime(&today);
    os_ << "\n%\n% Parameter set\n%\n";
    os_ << "L = " << param.L << "; % Area size\n";
    os_ << "v = " << param.v << "; % Speciation rate\n";
    os_ << "phi = " << param.phi << "; % Strength of phylogenetic effect\n";
    os_ << "psi = " << param.psi << "; % Strength of abundance effect\n";
    os_ << "sigmaA = " << param.sigmaA << "; % Abundance standard deviation\n";
    os_ << "sigmaB = " << param.sigmaB << "; % Dispersal standard deviation\n";
    os_ << "implicit = " << (param.implicit ? "true" : "false") << "; % Spacial implicit model\n";
    os_ << "dominant = " << (param.dominant ? "true" : "false") << "; % Dominant mode\n";
    os_ << "turnover = " << param.ticks << "; % Turnovers\n";
    os_ << "totalTicks = " << param.totalTicks << "; % total ticks\n";
    os_ << "log_interval = " << param.log_interval << "; % Log interval for M, D and R [events]\n";
    os_ << EventHeader;
    os_ << "events = [];\n";
    os_ << "extinctions = [];\n";
    os_ << "speciations = [];\n";
    os_ << "\n%\n% Snapshot log of M, D, and R\n";
    os_ << "% the first and the last set in a simulation are always logged\n";
    os_ << "% use M{end}, D{end} and R{end} to get the last records\n%\n";
    os_ << "M = {};    % cell array of M-matrices\n";
    os_ << "D = {};    % cell array of D-matrices\n";
    os_ << "R = {};    % cell array of abundance histograms\n";
    os_ << "sT = [];   % vector of snapshot times\n";
    os_ << "\n\n%!!!!!!!!!!!!! DO NOT EDIT AFTER THIS LINE !!!!!!!!!!!!!\n\n";
    os_ << "% " << cpp_skip_label_read_back << std::endl;
  }


  MatLogger::~MatLogger()
  {
  }


  void MatLogger::logEvents(std::vector<event> const& events)
  {
    os_ << "% " << cpp_skip_label_event_log << '\n';
    os_ << "\nevents = [\n";
    for (auto const& e : events) { os_ << e << '\n'; }
    os_ << "];\n";
    os_ << "if (length(events) > 0)\n";
    os_ << "    extinctions = events(events(:, 6) == -1, :);\n";
    os_ << "    speciations = events(events(:, 6) >= 0, :);\n";
    os_ << "end;" << std::endl;
  }


  void MatLogger::logSnapshot(int64_t T, square_buffer<int> const& M, GDM const& D, std::vector<int> const& R)
  {
    os_ << "M{length(M)+1} = [\n";
    os_ << M << "];" << std::endl;
    os_ << "D{length(D)+1} = [\n";
    os_ << D.data() << "];" << std::endl;
    os_ << "R{length(R)+1} = [ ";
    for (size_t i = 0; i < R.size(); ++i) { os_ << R[i] << ' '; }
    os_ << "];" << std::endl;
    os_ << "sT = [sT; " << T << "];" << std::endl;
  }


  void MatLogger::logLastSnapshot(int64_t T, square_buffer<int> const& M, GDM const& D, std::vector<int> const& R)
  {
    os_ << "% " << cpp_skip_label_last_log << '\n';
    logSnapshot(T, M, D, R);
  }
    

  void MatLogger::logState(Parameter const& param)
  {
    os_ << "elapsedTime = " << param.elapsed_time << ";\n\n";
    os_ << "\n%{\n";
    os_ << cpp_skip_label_state << '\n';
    os_ << '"' << JCVersion << "\"\n";
    os_ << param << '\n';
    os_ << "%}" << std::endl;
  }

} // namespace jc

