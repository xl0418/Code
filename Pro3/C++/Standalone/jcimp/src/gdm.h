/*! \file gdm.h
 *  \brief Genetic distance matrix
 *  \author Hanno Hildenbrandt
 */

#ifndef JC_GDM_H_INCLUDED
#define JC_GDM_H_INCLUDED

#include <cassert>
#include <cstdint>
#include "square_buffer.h"


namespace jc { 
  

  /// Genetic Distance Matrix
  class GDM
  {
  public:
    typedef square_buffer<int64_t> buffer;

  public:
    typedef buffer::value_type value_type;

    /// Creates the initial 1x1 zero-matrix
    GDM();

    /// Creates the GDM from square_buffer
    GDM(square_buffer<value_type> const& buf);

    void update() { ++dt_; }                  ///< Increase off-diagonal elements (tracked)
    void erase(int speciesId);                ///< Erase species \p speciesId from GDM
    void specification(int ancestorId);       ///< Handle specification
    
    value_type sum() const { return sum_ + dt_ * (buf_.elements() - buf_.n()); }  ///< Returns sum(Gij)

    value_type lastSum() const 
    { 
      assert(dt_ > 0);
      return sum_ + (dt_ - 1) * (buf_.elements() - buf_.n());
    }
    
    ///< element access
    value_type operator()(int i, int j) const 
    {
      assert(j != i); 
      return buf_(i, j) + dt_; 
    }
    
    value_type Tij(int i, int j) const 
    {
      assert(j != i); 
      assert(dt_ > 0);
      return buf_(i, j) + (dt_ - 1); 
    }
    
    buffer asMatrix() const;                  ///< Returns tidy matrix

  private:
    void swap(GDM& rhs);
    void symCopy(buffer& dst);
    void applyDt(buffer& buf);
    void recalcSum();

    buffer buf_, tmp_;        // ping-pong buffer
    value_type sum_;          // tracked sum( GDM_ij )
    value_type dt_;           // ticks since last change of size
  };


} // namespace jc

#endif
