// neutral in R model; psi == 0

#include "jc_base.h"


namespace jc {


  class ImplicitModel : public ModelBase
  {
  public:
    ImplicitModel(Parameter const& param) : ModelBase(param)
    {
      powPsi_.reserve(param.N);
      for (int i = 0; i <= param.N; ++i) {
        powPsi_.push_back(std::pow(double(i), param_.psi));
      }
      refresh_cdf();
    };

  protected:
    void do_run(int64_t T0);
    void refresh_cdf();
    void refresh_cdf_partial(int sp, int Rold);
    void refresh_pdfR();
    void refresh_pdfR_partial(int sp, int Rold);

  private:
    std::vector<double> pdfR_, pdfD_;
    rndutils::mutable_discrete_distribution<> cdf_;
    std::vector<double> powPsi_;
  };


  void ImplicitModel::do_run(int64_t T0)
  {
    auto specificationDist = std::bernoulli_distribution(param_.v);
    auto speciesDist = rndutils::mutable_discrete_distribution<>(R_.begin(), R_.end());
    for (int64_t T = T0; T < param_.ticks; ++T)
    {
      bool log = false;                            // assume uneventful time step
      D_.update();

      // Choose one individual to die
      auto victim = speciesDist(RndEng);
      auto Rold = R_[victim];
      if (Rold == 1)
      { // we have an extinction event
        log = true;
        D_.erase(victim);                          // remove victim species from GDM
        R_.erase(R_.begin() + victim);             // remove victim species from abundance vector
        speciesDist.mutate(R_);                    // update species sampler (full update)
        events_.push_back(event{T, R_.size(), victim, -1});
        refresh_cdf();
      }
      else
      {
        auto Rold = R_[victim]--;                  // adjust abundance
        speciesDist.mutate_partial(R_, victim);    // update species sampler (partial update)
        refresh_cdf_partial(victim, Rold);
      }

      // choose intruder
      auto sp = cdf_(RndEng);
      if (specificationDist(RndEng))
      { // we have an specification event
        log = true;
        auto ancestor = sp;
        D_.specification(sp);                      // insert sp in D
        R_.push_back(1);                           // insert sp in R
        sp = static_cast<int>(R_.size() - 1);      // sp : new species
        speciesDist.mutate(R_);                    // update species sample (full update)
        events_.push_back(event{T, R_.size(), sp, ancestor});
      }
      else
      {
        Rold = R_[sp]++;                           // adjust abundance
        speciesDist.mutate_partial(R_, sp);        // update species sampler (partial update)
      }
      log ? refresh_cdf() : refresh_cdf_partial(sp, Rold);

      handleOutput(log, T);
    }
  }


  void ImplicitModel::refresh_pdfR_partial(int sp, int Rold)
  {
    const int n = static_cast<int>(R_.size());
    // pdfR_[sp] must be completely recalculated
    double pA = -1.0;       // remove the i == sp summand in the loop below
    double powR = powPsi_[R_[sp]];
    for (int i = 0; i < n; ++i)
    {
      pA += powPsi_[R_[i]] / powR;
    }
    pdfR_[sp] = pA;
    // pdfR_[j], j != sp must be partially recalculated
    double powRold = powPsi_[Rold];
    for (int j = 0; j < sp; ++j)
    {
      auto powRj = powPsi_[R_[j]];
      double pA = powR / powRj;
      double pAold = powRold / powRj;
      pdfR_[j] += pA - pAold;
    }
    for (int j = sp+1; j < n; ++j)
    {
      auto powRj = powPsi_[R_[j]];
      double pA = powR / powRj;
      double pAold = powRold / powRj;
      pdfR_[j] += pA - pAold;
    }
  }
  
  
  void ImplicitModel::refresh_pdfR()
  {
    const int n = static_cast<int>(R_.size());
    pdfR_.assign(n, 0.0);
    for (int j = 0; j < n; ++j)
    {
      double powRj = powPsi_[R_[j]];
      for (int i = j + 1; i < n; ++i)
      {
        double pA = powPsi_[R_[i]] / powRj;
        pdfR_[j] += pA;
        pdfR_[i] += 1.0 / pA;
      }
    }
  }

  
  void ImplicitModel::refresh_cdf()
  {
    refresh_pdfR();
    const int n = static_cast<int>(R_.size());
    pdfD_.assign(n, 0.0);
    const double logDsum = std::log(D_.sum());
    const double phi = param_.phi;
    for (int j = 0; j < n; ++j)
    {
      for (int i = j + 1; i < n; ++i)
      {
        auto logDij = std::log(D_(i, j));   // Dij is guarantied to be non-zero
        auto x = logDij - logDsum;
        double pD = std::exp(phi * x);
        pdfD_[j] += pD;
        pdfD_[i] += pD;
      }
    }
    for (int i = 0; i < n; ++i) 
    {
      pdfD_[i] *= R_[i] * pdfR_[i];
    }
    cdf_.mutate(pdfD_);
  }


  void ImplicitModel::refresh_cdf_partial(int sp, int Rold)
  {
    refresh_pdfR_partial(sp, Rold);
    const int n = static_cast<int>(R_.size());
    const double phi = param_.phi;
    auto logDsum = std::log(D_.sum());
    pdfD_[sp] = 0;
    for (int i = 0; i < sp; ++i)
    {
      auto logDij = std::log(D_(i, sp));
      auto x = logDij - logDsum;
      double pD = std::exp(phi * x);
      pdfD_[sp] += R_[sp] * pD * pdfR_[sp];
    }
    for (int i = sp+1; i < n; ++i)
    {
      auto logDij = std::log(D_(i, sp));
      auto x = logDij - logDsum;
      double pD = std::exp(phi * x);
      pdfD_[sp] += R_[sp] * pD * pdfR_[sp];
    }
    cdf_.mutate(pdfD_);
  }


  ModelBase* CreateImplicitModel(Parameter const& param)
  {
    return new ImplicitModel(param);
  }

}
