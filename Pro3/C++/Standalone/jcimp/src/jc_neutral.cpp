// neutral model: phi == 0m psi == 0

#include "jc_base.h"


namespace jc {

  class NeutralModel : public ModelBase
  {
  public:
    NeutralModel(Parameter const& param) : ModelBase(param)
    {
      assert(param.phi == 0 && param.psi == 0.0 && "NeutrralModel: phi and psi shall be zero");
    };

  protected:
    void do_run(int64_t T0);

  private:
  };


  void NeutralModel::do_run(int64_t T0)
  {
    auto specificationDist = std::bernoulli_distribution(param_.v);
    auto speciesDist = rndutils::mutable_discrete_distribution<>(R_.begin(), R_.end());
    for (int64_t T = T0; T < param_.ticks; ++T)
    {
      bool log = false;                            // assume uneventful time step
      D_.update();

      // Choose one individual to die
      auto victim = speciesDist(RndEng);
      if (R_[victim] == 1)
      { // we have an extinction event
        log = true;
        D_.erase(victim);                             // Remove species victim from GDM
        R_.erase(R_.begin() + victim);                // Remove abundance entry for species victim
        speciesDist.mutate(R_);
        events_.push_back(event{T, R_.size(), victim, -1});
      }
      else
      {
        --R_[victim];
        speciesDist.mutate_partial(R_, victim);
      }

      // choose intruder
      auto sp = speciesDist(RndEng);
      if (specificationDist(RndEng))
      { // we have an specification event
        log = true;
        auto ancestor = sp;
        D_.specification(sp);
        R_.push_back(1);
        sp = static_cast<int>(R_.size() - 1);
        speciesDist.mutate(R_);
        events_.push_back(event{T, R_.size(), sp, ancestor});
      }
      else
      {
        ++R_[sp];
        speciesDist.mutate_partial(R_, sp);
      }

      handleOutput(log, T);
    }
  }


  ModelBase* CreateNeutralModel(Parameter const& param)
  {
    return new NeutralModel(param);
  }

}
