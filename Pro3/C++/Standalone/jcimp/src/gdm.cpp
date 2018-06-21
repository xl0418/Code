#include <cassert>
#include <utility>
#include <numeric>
#include <algorithm>
#include "gdm.h"


namespace jc {
  
  using detail::blk_copy;


  GDM::GDM() 
    : buf_(1)
  {
    sum_ = 1;
    dt_ = 0;
  }


  GDM::GDM(square_buffer<value_type> const& buf) : buf_(buf)
  {
    recalcSum();
    dt_ = 0;
  }


  void GDM::erase(int speciesId)
  {
    applyDt(buf_); dt_ = 0;
    // create a copy without the row & column speciesId
    const int n = static_cast<int>(buf_.n());
    tmp_.resize(buf_.n() - 1);
    int r = n - speciesId - 1;
    blk_copy(tmp_, buf_, speciesId, speciesId, 0, 0, 0, 0);                           // left top
    blk_copy(tmp_, buf_, r, speciesId, speciesId, 0, speciesId + 1, 0);               // right top
    blk_copy(tmp_, buf_, speciesId, r, 0, speciesId,  0, speciesId + 1);              // left bottom
    blk_copy(tmp_, buf_, r, r, speciesId, speciesId, speciesId + 1, speciesId + 1);   // right bottom
    buf_.swap(tmp_);
    recalcSum();
  }


  void GDM::specification(int ancestorId)
  {
    applyDt(buf_); dt_ = 0;
    const int n = static_cast<int>(buf_.n());
    tmp_.resize(buf_.n() + 1);
    symCopy(tmp_);                                              // copy matrix
    blk_copy(tmp_, tmp_, n, 1, 0, n, 0, ancestorId);            // duplicate row into last row
    blk_copy(tmp_, tmp_, 1, n, n, 0, ancestorId, 0);            // duplicate column into last column
    tmp_(n, n) = 0;                                             // distance to itself <- 0
    tmp_(n, ancestorId) = tmp_(ancestorId, n) = 1;              // distance to ancestor <- 1
    buf_.swap(tmp_);
    recalcSum();
  }


  // copy buf_[upperTri] -> tmp[upperTri], tmp[lowerTri]
  void GDM::symCopy(buffer& dst)
  {
    auto const n = static_cast<int>(buf_.n());
    auto const nt = static_cast<int>(dst.n());
    const GDM::buffer::value_type* __restrict p = &(*buf_.begin());
    GDM::buffer::value_type* __restrict pt = &(*dst.begin());
    for (int i = 0; i < n; ++i) {
      for (int j= i + 1; j < n; ++j) {
        pt[i * nt + j] = pt[i + nt * j] = p[i * n + j];
      }
    }
  }


  void GDM::applyDt(GDM::buffer& buf)
  {
    GDM::buffer::value_type* __restrict p = &(*buf.begin());
    auto const n = static_cast<int>(buf.n());
    for (int i = 0; i < n - 1; ++i) {
      for (int j= i + 1; j < n; ++j) {
        p[i * n + j] += dt_;
        p[i + n * j] += dt_;
      }
    }
  }
  

  void GDM::recalcSum()
  {
    sum_ = 0;
    const GDM::buffer::value_type* __restrict p = &(*buf_.cbegin());
    auto const n = static_cast<int>(buf_.n());
    for (int i = 0; i < n - 1; ++i) {
      for (int j= i + 1; j < n; ++j) {
        sum_ += p[i * n + j];
      }
    }
    sum_ *= 2;
  }
  

  GDM::buffer GDM::asMatrix() const
  {
    buffer tmp(buf_.n());
    const_cast<GDM*>(this)->symCopy(tmp);
    const_cast<GDM*>(this)->applyDt(tmp);
    return tmp;
  }
  

  void GDM::swap(GDM& rhs)
  {
    buf_.swap(rhs.buf_);
    tmp_.swap(rhs.tmp_);
    std::swap(sum_, rhs.sum_);
  }


} // namespace jc
