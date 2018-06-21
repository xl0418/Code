#include <stdexcept>
#include "rndutils.hpp"
#include <ctime>
#include <string>
#include "jc_base.h"
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
% events_(:,1)    The time the event occurred.
% events_(:,2)    Number of species after the event.
% events_(:,3:4)  Position on the grid.
% events_(:,5)    New or extinct species.
% events_(:,6)    Ancestor in case of speciation event, -1 otherwise.
%
% extinctions    Extinction events_ from events_
% speciations    Speciation events_ from events_
%
)m";


  std::ostream& operator<<(std::ostream& os, event const& e)
  {
    os << "  " << e.T << ' ' << e.NS << ' ' << e.sp << ' ' << e.ancestor << ';';
    return os;
  }


  std::ostream& operator<<(std::ostream& os, Parameter const& p)
  {
    os << p.N << ' ' << p.v << ' ';
    os << p.phi << ' ' << p.psi << ' ';
    os << p.ticks << ' ';
    os << '"' << p.filename.c_str() << "\" ";
    os << p.log_interval;
    return os;
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
    os_ << "N = " << param.N << "; % Individuals\n";
    os_ << "v = " << param.v << "; % Speciation rate\n";
    os_ << "phi = " << param.phi << "; % Strength of phylogenetic effect\n";
    os_ << "psi = " << param.psi << "; % Strength of abundance effect\n";
    os_ << "implicit = " << (param.implicit ? "true" : "false") << "; % Spacial implicit model\n";
    os_ << "turnover = " << param.ticks << "; % Turnovers\n";
    os_ << "log_interval = " << param.log_interval << "; % Log interval for M, D and R [events_]\n";
    os_ << "repetition = " << param.rep << "; % Repetition\n";

    os_ << EventHeader;
    os_ << "events = [];\n";
    os_ << "extinctions = [];\n";
    os_ << "speciations = [];\n";
    os_ << "\n%\n% Snapshot log of M, D, and R\n";
    os_ << "% the first and the last set in a simulation are always logged\n";
    os_ << "% use D{end} and R{end} to get the last records\n%\n";
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
    os_ << "    extinctions = events(events(:, 4) == -1, :);\n";
    os_ << "    speciations = events(events(:, 4) >= 0, :);\n";
    os_ << "end;" << std::endl;
  }


  void MatLogger::logSnapshot(int64_t T, GDM const& D, std::vector<int> const& R)
  {
    os_ << "D{length(D)+1} = [\n";
    os_ << D.asMatrix() << "];" << std::endl;
    os_ << "R{length(R)+1} = [ ";
    for (size_t i = 0; i < R.size(); ++i) { os_ << R[i] << ' '; }
    os_ << "];" << std::endl;
    os_ << "sT = [sT; " << T << "];" << std::endl;
  }


  void MatLogger::logLastSnapshot(int64_t T, GDM const& D, std::vector<int> const& R)
  {
    os_ << "% " << cpp_skip_label_last_log << '\n';
    logSnapshot(T, D, R);
  }
    

  void MatLogger::logState(Parameter const& param)
  {
    os_ << "\n%{\n";
    os_ << cpp_skip_label_state << '\n';
    os_ << '"' << JCVersion << "\"\n";
    os_ << RndEng << '\n';
    os_ << param << '\n';
    os_ << "%}" << std::endl;
  }

} // namespace jc

