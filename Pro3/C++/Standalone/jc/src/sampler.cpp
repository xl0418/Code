#include <iostream>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <algorithm>
#include "rndutils.hpp"
#include <cassert>
#include "sampler.h"


namespace jc {


  class ExplicitSampler : public Sampler
  {
  public:
    ExplicitSampler(Parameter const& param);
    ~ExplicitSampler() {}

    int apply(std::pair<int, int> const& pos,
              square_buffer<int> const& M,
              GDM const& D,
              std::vector<int> const& R,
              double rnd);

  private:
    const double phi_;
    const double psi_;
    const double factA_;
    const double factB_;
    std::vector<int> dd_;
    std::vector<int> scan_;
    std::vector<expr::expr> exprRa_;
    std::vector<expr::expr> exprRb_;
    std::vector<expr::expr> exprP_;
    expr::memorizer ma_;
    expr::memorizer mb_;
    std::vector<int> minDD_;
    std::vector<double> cdf_;
  };


  ExplicitSampler::ExplicitSampler(Parameter const& param)
  : phi_(param.phi),
    psi_(param.psi),
    factA_(-0.5 / (param.sigmaA * param.sigmaA)),
    factB_(-0.5 / (param.sigmaB * param.sigmaB)),
    dd_(param.L * param.L + 1),
    ma_(param.L, factA_),
    mb_(param.L, factB_)
  {
  }


  int ExplicitSampler::apply(std::pair<int, int> const& pos,
                     square_buffer<int> const& M,
                     GDM const& D,
                     std::vector<int> const& R,
                     double rnd)
  {
    const int s = pos.first;
    const int t = pos.second;
    const size_t NS = R.size();         // number of species
    exprRa_.resize(NS);
    exprRb_.resize(NS);
    minDD_.assign(NS, std::numeric_limits<int>::max());

    // setup of the distance-scan
    scan_.resize(NS + 1);
    int cs = scan_[0] = 0;
    for (size_t i = 0; i < NS; ++i)
    {
      cs += R[i];
      scan_[i + 1] = cs;
    }

    // scan squared distances
    {
      const int L = static_cast<int>(M.n());
      int left = std::max(0, s - L);
      int right = std::min(L - 1, s + L);
      int top = std::max(0, t - L);
      int bottom = std::min(L - 1, t + L);
      for (int y = top; y <= bottom; ++y)
      {
        int ddy = (y - t) * (y - t);
        for (int x = left; x <= right; ++x)
        {
          int dd = ddy + (s - x) * (s - x);
          if (dd)
          {
            int sp = M(x, y);
            dd_[scan_[sp]++] = dd;
            if (minDD_[sp] > dd) minDD_[sp] = dd;
          }
        }
      }
    }
    size_t first = 0;
    for (size_t i = 0; i < NS; ++i)
    {
      auto tmp = expr::powsum(&dd_[first], &dd_[scan_[i]], ma_, mb_, minDD_[i]);
      exprRa_[i] = tmp.first;
      exprRb_[i] = tmp.second;
      first += R[i];
    }
    expr::normalize(exprRa_.begin(), exprRa_.end());
    expr::normalize(exprRb_.begin(), exprRb_.end());

    cdf_.resize(NS);
    auto const& Dbuf = D.data();
    const auto Dsum = D.sum();
    for (size_t j = 0; j < NS; ++j)
    {
      exprP_.clear();
      for (size_t i = 0; i < NS; ++i)
      {
        // A bit of caution in the calculation of pow(D, phi):
        // The elements of D could become very large, and therefore
        // divisions by their sum might underflow.
        if (j != i)
        {
          auto Dij = Dbuf(i, j);              // guarantied to be non-zero
          // calculate x^phi as (1/x)^-phi
          auto inorm = std::div(Dsum, Dij);   // 1/norm
          double pD = std::pow(inorm.quot + static_cast<double>(inorm.rem) / static_cast<double>(Dij), -phi_);
          exprP_.push_back(pD * expr::pow(exprRa_[i] / exprRa_[j], psi_) * exprRb_[j]);
        }
        else
        {
          exprP_.push_back(exprRb_[j]);
        }
      }
      cdf_[j] = R[j] * expr::sum(exprP_.cbegin(), exprP_.cend()).value();
    }
    double val = 0.0;
    for (auto& x : cdf_) { val += x; x = val; }
    double threshold = rnd * cdf_.back();
    return static_cast<int>(std::count_if(cdf_.cbegin(), cdf_.cend(), [threshold](double x){ return x < threshold; }));
  }


  class ImplicitSampler : public Sampler
  {
  public:
    ImplicitSampler(Parameter const& param);
    ~ImplicitSampler() {}

    int apply(std::pair<int, int> const&,
              square_buffer<int> const&,
              GDM const& D,
              std::vector<int> const& R,
              double rnd);

  private:
    const double phi_;
    const double psi_;
    std::vector<double> cdf_;
    square_buffer<double> P_;
  };


  ImplicitSampler::ImplicitSampler(Parameter const& param)
  : phi_(param.phi),
    psi_(param.psi)
  {
  }


  inline double fastPow(double a, double b) 
  {
    union {
      double d;
      int x[2];
    } u = { a };
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
  }

  int ImplicitSampler::apply(std::pair<int, int> const&,
                             square_buffer<int> const&,
                             GDM const& D,
                             std::vector<int> const& R,
                             double rnd)
  {
    const size_t NS = R.size();         // number of species
    cdf_.resize(NS);
    auto const& Dbuf = D.data();
    const auto log2Dsum = std::log2(D.sum());
    {
      P_.resize(NS);
      P_.assign(0.0);
      for (size_t j = 0; j < NS; ++j)
      {
        for (size_t i = j+1; i < NS; ++i)
        {
          // A bit of caution in the calculation of pow(D, phi):
          // The elements of D could become very large, and therefore
          // divisions by their sum might underflow. 
          auto log2Dij = std::log2(Dbuf(i, j));              // guarantied to be non-zero
          auto x = log2Dij - log2Dsum;
//          // calculate x^phi as (1/x)^-phi
//          auto inorm = std::div(Dsum, Dij);   // 1/norm
//          auto x = inorm.quot + static_cast<double>(inorm.rem) / static_cast<double>(Dij);
          double pD = std::exp2(phi_ * x);
          double pA = std::pow(static_cast<double>(R[i]) / R[j], psi_);
          P_(j, i) += pD * pA;
          P_(i, j) += pD / pA;    // pA_ij = 1 / pA_ji
        }
      }
      for (size_t j = 0; j < NS; ++j)
      {
        double P = 0.0;
        for (size_t i = 0; i < NS; ++i)
        {
          P += P_(j, i);
        }
        cdf_[j] = R[j] * P;
      }
    }
    double val = 0.0;
    for (auto& x : cdf_) { val += x; x = val; }
    double threshold = rnd * cdf_.back();
    auto it = std::lower_bound(cdf_.cbegin(), cdf_.cend(), threshold);
    return static_cast<int>(std::distance(cdf_.cbegin(), it));
  }


  class ImplicitNeutralSampler : public Sampler
  {
  public:
    ImplicitNeutralSampler(Parameter const& param) {}
    ~ImplicitNeutralSampler() {}

    int apply(std::pair<int, int> const&,
              square_buffer<int> const&,
              GDM const& D,
              std::vector<int> const& R,
              double rnd);
  private:
    std::vector<double> cdf_;
  };


  int ImplicitNeutralSampler::apply(std::pair<int, int> const&,
                                    square_buffer<int> const&,
                                    GDM const&,
                                    std::vector<int> const& R,
                                    double rnd)
  {
    cdf_.resize(R.size());
    double val = 0.0;
    for (size_t i=0; i < R.size(); ++i) 
    {
      cdf_[i] = val += static_cast<double>(R[i]);
    }
    double threshold = rnd * cdf_.back();
    auto it = std::lower_bound(cdf_.cbegin(), cdf_.cend(), threshold);
    return static_cast<int>(std::distance(cdf_.cbegin(), it));
  }


  Sampler* CreateExplicitSampler(Parameter const& param)
  {
    return new ExplicitSampler(param);
  }


  Sampler* CreateImplicitSampler(Parameter const& param)
  {
    if ((param.phi == 0.0) && (param.psi == 0.0))
    {
      return new ImplicitNeutralSampler(param);
    }
    return new ImplicitSampler(param);
  }


  Sampler* Sampler::create(Parameter const& param)
  {
    if (param.implicit)
    {
      return CreateImplicitSampler(param);
    }
    else
    {
      return CreateExplicitSampler(param);
    }
  }


} // namespace jc
