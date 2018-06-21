#ifndef JC_BASE_H_INCLUDED
#define JC_BESE_H_INCLUDED


#include <iostream>
#include "jc.h"
#include "gdm.h"
#include "iomatlab.h"
#include "neutral.h"


namespace jc {


  extern rndutils::xorshift128 thread_local RndEng;
  

  class ModelBase
  {
    ModelBase(ModelBase const&) = delete;
    ModelBase& operator=(ModelBase const&) = delete;

  public:
    ModelBase(Parameter const& param)
    : param_(param),
      R_(jc::R0),
      D_(jc::D0),
      logger_(param)
    {
    }

    void run();

  protected:
    void handleOutput(bool log, int64_t T);
    virtual void do_run(int64_t T0) = 0;

    Parameter param_;
    std::vector<int> R_;                             // Abundance vector
    GDM D_;                                          // Genetic distance matrix    
    std::vector<event> events_;                      // Event log

  private:
    MatLogger logger_;                               // Matlab output
  };


  inline void ModelBase::run()
  {
    if (!param_.implicit) throw std::runtime_error("This is the implicit JC model. Run with the '--implicit' option");
    auto tic = std::chrono::high_resolution_clock::now();
    do_run(0);
    handleOutput(false, param_.ticks);
    auto toc = std::chrono::high_resolution_clock::now();
    auto dt = double(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count());
    if (param_.verbose) std::cout << dt << " ms  " << 1000.0 * dt / param_.ticks << " us per turnover" << std::endl;
    // Epilogue
    logger_.logLastSnapshot(param_.ticks, D_, R_);
    logger_.logEvents(events_);
    logger_.logState(param_);
  }


  inline void ModelBase::handleOutput(bool log, int64_t T)
  {
    if (log && (T > param_.ticks >> 1) && (events_.size() % param_.log_interval == 0))
    {
      logger_.logSnapshot(T, D_, R_);
    }
    if (T % 1000000 == 0 && param_.verbose)
    {
      std::cout << '[' << param_.rep << "] " << param_.phi << ' ' << param_.psi << "  " << T << "  " << R_.size();
      for (size_t i = 0; i < 15 && i < R_.size(); ++i) std::cout << ' ' << R_[i];
      std::cout << '\n';
    }
  }
   
}


#endif

