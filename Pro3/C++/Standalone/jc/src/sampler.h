/*! \file sampler.h
*  \brief Handles Gaussian kernel
*  \author Hanno Hildenbrandt
*/

#ifndef JANZENCONNELL_SAMPLER_H_INCLUDED
#define JANZENCONNELL_SAMPLER_H_INCLUDED

#include <vector>
#include <utility>
#include "jc.h"
#include "gdm.h"
#include "exp_expr.h"


namespace jc { 
  

  class Sampler
  {
    Sampler(Sampler const&) = delete;
    Sampler& operator=(Sampler const&) = delete;

  protected:
    Sampler() {}

  public:
    static Sampler* create(Parameter const& param);
    virtual ~Sampler() {}
    virtual int apply(std::pair<int, int> const& pos,
                      square_buffer<int> const& M,
                      GDM const& D, 
                      std::vector<int> const& R,
                      double rnd) = 0;
  };


} // namespace jc

#endif
