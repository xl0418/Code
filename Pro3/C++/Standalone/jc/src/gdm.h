/*! \file gdm.h
 *  \brief Genetic distance matrix
 *  \author Hanno Hildenbrandt
 */

#ifndef JC_GDM_H_INCLUDED
#define JC_GDM_H_INCLUDED

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

    /// Creates the initial 2x2 matrix.
    /// 
    /// Initialize the GDM to:
    ///   0 1
    ///   1 0
    GDM();

    /// Creates the GDM from square_buffer
    GDM(square_buffer<value_type> const& buf);

    void erase(int speciesId);                ///< Erase species \p speciesId from GDM
    void specification(int ancestorId);       ///< Handle specification
    void update(int64_t updateT);             ///< Increase off-diagonal elements by updateT
    value_type sum() const { return sum_; }   ///< Returns sum( GDM_ij )

    /// const access as linear array
    buffer const& data() const { return buf_; }

  private:
    void swap(GDM& rhs);

    buffer buf_, tmp_;    // ping-pong buffer
    value_type sum_;      // tracked sum( GDM_ij )
  };


} // namespace jc

#endif
