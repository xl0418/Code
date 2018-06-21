/*! \file square_buffer.h
*  \brief 2D memory buffer
*  \author Hanno Hildenbrandt
*/

#ifndef JC_SQUARE_BUFFER_H_INCLUDED
#define JC_SQUARE_BUFFER_H_INCLUDED

#include <vector>
#include <type_traits>


namespace jc {


  template <typename T>
  class square_buffer
  {
  public:
    typedef std::vector<T>                          linear_buffer;
    typedef typename linear_buffer::size_type       size_type;
    typedef typename linear_buffer::value_type      value_type;
    typedef typename linear_buffer::iterator        iterator;
    typedef typename linear_buffer::const_iterator  const_iterator;

    square_buffer() : n_(0) {}
    square_buffer(size_type n) : n_(n), buf_(n * n, T()) {}
    square_buffer(square_buffer const& buf) : n_(buf.n()), buf_(buf.buf_) {}

    void resize(size_type n) { buf_.resize(n * n, T()); n_ = n; }
    void assign(value_type val) { buf_.assign(buf_.size(), val);  }

    size_type n() const { return n_; }
    size_type elements() const { return n_ * n_; }

    // Access as linear array
    iterator begin() { return buf_.begin(); }
    const_iterator begin() const { return buf_.begin(); }
    const_iterator cbegin() const { return buf_.begin(); }
    iterator end() { return buf_.end(); }
    const_iterator end() const { return buf_.end(); }
    const_iterator cend() const { return buf_.end(); }

    // Element access
    value_type& operator()(size_type j) { return buf_[j]; }
    value_type operator()(size_type j) const { return buf_[j]; }
    value_type& operator()(size_type s,size_type t) { return buf_[t * n_ + s]; }
    value_type operator()(size_type s,size_type t) const { return buf_[t * n_ + s]; }

    void swap(square_buffer& other) 
    { 
      std::swap(n_, other.n_);
      buf_.swap(other.buf_); 
    }

  private:
    size_type n_;
    linear_buffer buf_;
  };


  namespace detail {


    template <typename T>
    inline void blk_copy(T* dst, T const* src, int w, int h, size_t dstS, size_t dstT, size_t srcS, size_t srcT, size_t dstStride, size_t srcStride)
    {
      if (w <= 0 || h <= 0) return;
      dst += dstT * dstStride + dstS;
      src += srcT * srcStride + srcS;
      for (int j = 0; j < h; ++j, dst += dstStride, src += srcStride)
      {
        for (int i = 0; i < w; ++i)
        {
          *(dst + i) = *(src + i);
        }
      }
    }


    template <typename T>
    inline void blk_copy(square_buffer<T>& dst, square_buffer<T> const& src, int w, int h, int dstS, int dstT, int srcS, int srcT)
    {
      blk_copy(&*dst.begin(), &*src.cbegin(), w, h, dstS, dstT, srcS, srcT, dst.n(), src.n());
    }


  } // namespace detail

} // namespace jc

#endif
