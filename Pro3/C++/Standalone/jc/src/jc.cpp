#include <iostream>
#include <memory>
#include <fstream>
#include "rndutils.hpp"
#include <numeric>
#include <algorithm>
#include <chrono>
#include "jc.h"
#include "sampler.h"
#include "gdm.h"
#include "iomatlab.h"


namespace jc {


  void print_event(event const& e, std::vector<int> const& R)
  {
    std::cout << e << "   ";
    for (size_t i = 0; i < std::min<size_t>(10u, R.size()); ++i)
    {
      std::cout << R[i] << ' ';
    }
    std::cout << std::endl;
  }


  class Model
  {
    Model(Model const&) = delete;
    Model& operator=(Model const&) = delete;

  public:
    Model();
    Model(Parameter const& param);
    void run();

  private:
    bool checkExitCondition(int64_t T);
    std::pair<int, int> randomPos();

    const Parameter param_;
    int64_t T_;                                      // Time step
    int64_t Td_;                                     // The time step equilibrium is reached
    rndutils::default_engine reng_;                              // Random number generator
    std::uniform_real_distribution<double> rnd01_;   // Uniform distribution [0, 1)
    std::uniform_int_distribution<int> rndJ_;        // Uniform distribution [0, J)
    square_buffer<int> M_;                           // Grid
    std::vector<int> R_;                             // Abundance vector
    GDM D_;                                          // Genetic distance matrix    
    std::vector<event> events;                       // Event log
  };


  Model::Model(Parameter const& param)
    : param_(param),
    T_(0), Td_(0),
    reng_(rndutils::make_random_engine<>()),
    rnd01_(0.0, 1.0),
    rndJ_(0, param.L * param.L - 1),      // L*L-1 : uniform_int_distribution is inclusive
    M_(param.L),
    R_({ param.L * param.L - 1, 1 })
  {
    // place singleton of species 1
    auto pos = randomPos();
    M_(pos.first, pos.second) = 1;
  }


  bool Model::checkExitCondition(int64_t T)
  {
    if (param_.dominant)
    {
      if (Td_ == 0) 
      {
        if (R_[0] != *std::max(R_.begin(), R_.end()))
        {
          // we are approaching equilibrium
          Td_ = T;
        }
      }
      return (Td_ == 0) ? true : (T < (Td_ + param_.ticks));
    }
    return T < param_.ticks;
  }


  void Model::run()
  {
    auto tic = std::chrono::system_clock::now();

    MatLogger logger(param_);
    std::unique_ptr<Sampler> sampler(Sampler::create(param_));
    auto ptic = std::chrono::high_resolution_clock::now();        // profiler clock
    int64_t updateT = 0;
    for (int64_t T = T_; checkExitCondition(T); ++T, ++param_.totalTicks)
    {
      bool log = false;                            // assume uneventful time step

      ++updateT;

      // Choose one individual to die
      auto pos = randomPos();
      int sp = M_(pos.first, pos.second);

      if (R_[sp] == 1)
      {
        // we have an extinction event
        D_.update(updateT);
        updateT = 0;
        log = true;
        D_.erase(sp);                             // Remove species sp from GDM
        R_.erase(R_.begin() + sp);                // Remove abundance entry for species sp
        for (auto& x : M_) { if (x > sp) --x; }   // Adjust species id
        events.push_back(event{T, R_.size(), pos, sp, -1});
        if (param_.verbose) print_event(events.back(), R_);
      }
      else
      {
        --R_[sp];
      }

      // choose intruder species
      sp = sampler->apply(pos, M_, D_, R_, rnd01_(reng_));

      if (param_.v > rnd01_(reng_))
      {
        // we have an specification event
        if (!log)
        { // D-matrix not yet updated
          D_.update(updateT);
          updateT = 0;
          log = true;
        }
        auto ancestor = sp;
        D_.specification(sp);
        R_.push_back(1);
        sp = static_cast<int>(R_.size()) - 1;
        events.push_back(event{T, R_.size(), pos, sp, ancestor});
        if (param_.verbose) print_event(events.back(), R_);
      }
      else
      {
        ++R_[sp];
      }
      M_(pos.first, pos.second) = sp;

      // Log snapshot if necessary
      if ((T == T_) || (T < param_.log_first) || (log && ((param_.log_interval > 0) && (events.size() % param_.log_interval == 0))))
      {
        logger.logSnapshot(T, M_, D_, R_);
      }
      if (T % 1000000 == 0)
      {
        std::cout << T << ' ' << R_.size() << ' ' << R_[0] << ' ' << R_[1] << ' ' << D_.sum() << '\n';
      }
    }
    auto ptoc = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(ptoc - ptic);
    std::cout << milli.count() << ' ' << milli.count() / double(param_.ticks) << " ms" << std::endl;

    // Epilogue
    auto toc = std::chrono::system_clock::now();
    logger.logLastSnapshot(T_, M_, D_, R_);
    logger.logEvents(events);
    std::chrono::duration<double> elapsed = toc - tic;
    const_cast<Parameter&>(param_).elapsed_time += elapsed.count();
    logger.logState(param_);
  }


  inline std::pair<int, int> Model::randomPos()
  {
    auto const linPos = rndJ_(reng_);
    return{ linPos / param_.L, linPos % param_.L };
  }


  void RunModel(Parameter const& param)
  {
    Model model(param);
    model.run();
  }


} // namespace jc
