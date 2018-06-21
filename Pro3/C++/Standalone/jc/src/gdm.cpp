#include <utility>
#include <numeric>
#include <algorithm>
#include "gdm.h"


namespace jc {
  
  using detail::blk_copy;


  GDM::GDM() 
    : buf_(2)
  {
    update(1);
    sum_ = 2;
  }


  GDM::GDM(square_buffer<value_type> const& buf) : buf_(buf)
  {
    sum_ = std::accumulate(buf_.cbegin(), buf_.cend(), value_type(0));      // calculate sum
  }


  void GDM::erase(int speciesId)
  {
    // create a copy without the row & column speciesId
    const int n = static_cast<int>(buf_.n());
    tmp_.resize(buf_.n() - 1);
    int r = n - speciesId - 1;
    blk_copy(tmp_, buf_, speciesId, speciesId, 0, 0, 0, 0);                           // left top
    blk_copy(tmp_, buf_, r, speciesId, speciesId, 0, speciesId + 1, 0);               // right top
    blk_copy(tmp_, buf_, speciesId, r, 0, speciesId,  0, speciesId + 1);              // left bottom
    blk_copy(tmp_, buf_, r, r, speciesId, speciesId, speciesId + 1, speciesId + 1);   // right bottom
    buf_.swap(tmp_);
    sum_ = std::accumulate(buf_.cbegin(), buf_.cend(), value_type(0));      // recalc sum
  }


  void GDM::specification(int ancestorId)
  {
    const int n = static_cast<int>(buf_.n());
    tmp_.resize(buf_.n() + 1);
    blk_copy(tmp_, buf_, n, n, 0, 0, 0, 0);                     // copy matrix
    blk_copy(tmp_, buf_, n, 1, 0, n, 0, ancestorId);            // duplicate row into last row
    blk_copy(tmp_, buf_, 1, n, n, 0, ancestorId, 0);            // duplicate column into last column
    tmp_(n, n) = 0;                                             // distance to itself <- 0
    tmp_(n, ancestorId) = tmp_(ancestorId, n) = 1;              // distance to ancestor <- 1
    buf_.swap(tmp_);
    sum_ = std::accumulate(buf_.cbegin(), buf_.cend(), value_type(0));    // recalc sum
  }


  void GDM::update(int64_t updateT)
  {
    for (auto& x : buf_) { x += updateT; };
    const size_t n = buf_.n();
    for (size_t d = 0; d < n; ++d)
    {
      buf_(d, d) = 0;
    }
    // keep track of sum
    sum_ += value_type(buf_.elements() - buf_.n());
  }
  

  void GDM::swap(GDM& rhs)
  {
    buf_.swap(rhs.buf_);
    tmp_.swap(rhs.tmp_);
    std::swap(sum_, rhs.sum_);
  }


} // namespace jc
