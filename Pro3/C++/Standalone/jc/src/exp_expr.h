/*! \file exp_expr.h
*  \brief deferred exponential expressions assuming ridiculous small exponents
*  \author Hanno Hildenbrandt
*/

#ifndef JC_EXP_EXPR_H_INCLUDED
#define JC_EXP_EXPR_H_INCLUDED

#include <utility>
#include <numeric>
#include <algorithm>
#include <vector>


namespace expr {


  /// Holder of precomputed exp(y * dd) factors
  class memorizer
  {
  public:
    memorizer() {}

    /// Creates memorizer
    ///
    /// \param L grid size
    /// \param y exponent
    ///
    /// claculates exp(y * dd) for all posible \p dd = distance^2 on
    /// a L x L grid.
    ///
    memorizer(int L, double y);

    double y() const { return y_; }                  ///< returns common factor \p y in exp(y * dd)
    double val(int dd) const { return buf_[dd]; }    ///< returns exp(y * dd)

  private:
    std::vector<double> buf_;
    double y_;
  };


  /// \brief expression class for exp(n) * fact
  struct expr
  {
    double n;
    double fact;

    double value() const { return std::exp(n) * fact; }
    expr& operator *= (double y) { fact *= y; return *this; }
    expr& operator *= (expr const& y) { n += y.n; fact *= y.fact; return *this; }
    expr& operator /= (double y) { fact /= y; return *this; }
    expr& operator /= (expr const& y) { n -= y.n; fact /= y.fact; return *this; }
  };


  inline expr operator * (double y, expr const& x) { return{ x.n, x.fact * y }; }
  inline expr operator * (expr const& x, double y) { return{ x.n, x.fact * y }; }
  inline expr operator * (expr const& x, expr const& y) { return{ x.n + y.n, x.fact * y.fact }; }
  inline expr operator / (double y, expr const& x) { return{ x.n, x.fact / y }; }
  inline expr operator / (expr const& x, double y) { return{ x.n, x.fact / y }; }
  inline expr operator / (expr const& x, expr const& y) { return{ x.n - y.n, x.fact / y.fact }; }


  // x^a
  inline expr pow(expr const& x, double a)
  {
    return{ a * x.n, x.fact };
  }


  // (1/x)^a
  inline expr invpow(expr const& x, double a)
  {
    return{ -a * x.n, 1.0 / x.fact };
  }


  /// \deprecated left for debugging
  /// \brief Sums of exponential expression
  ///
  /// res.first <- exp(a * *first) + exp(a * *(first + 1)) + ... + exp(a * *(last -1))
  /// res.second <- exp(b * *first) + exp(b * *(first + 1)) + ... + exp(b * *(last -1))
  ///
  std::pair<expr, expr> powsum(const int* first, const int* last, double a, double b);


  /// \brief Sums of exponential expression
  ///
  /// \param minn minimal value in [first,last)
  ///
  /// same as \sa powsum above but uses information from the memorizes \p ma and \p mb
  ///
  std::pair<expr, expr> powsum(const int* first, const int* last, memorizer const& ma, class memorizer const& mb, int minn);


  /// sum of exp(n)*fact summands
  template <typename IT>
  inline expr sum(IT first, IT last)
  {
    double sum = 0.0;
    double maxn = std::max_element(first, last, [](expr const& a, expr const& b){ return a.n < b.n; })->n;
    for (; first != last; ++first)
    {
      sum += std::exp(first->n - maxn) * first->fact;
    }
    return{ maxn, sum };
  }


  /// normalize exp(n)*fact series
  template <typename IT>
  inline void normalize(IT first, IT last)
  {
    auto s = sum(first, last);
    for (; first != last; ++first)
    {
      *first /= s;
    }
  }

}

#endif
