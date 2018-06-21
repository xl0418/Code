/* exp_expr.cpp
 *
 * Copyright(c) 2015 Hanno Hildenbrandt
 *
 * This file is part of POJCE.
 *
 * POJCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include "exp_expr.h"


namespace expr {


  memorizer::memorizer(int L, double y) : y_(y)
  {
    for (int dd = 0; dd <= 2*L*L; ++dd)
    {
      buf_.push_back( std::exp(y_ * dd) );
    }
  }


  /// exp(a * *first) + exp(a * *(first + 1) + ... + exp(a * *(last-1)
  expr powsum(const int* first, const int* last, double a)
  {
    int minn = *std::min_element(first, last);
    double fact = 0.0;
    for (; first != last; ++first)
    {
      fact += std::exp(a * (*first - minn));
    }
    return{ a * minn, fact };
  }


  std::pair<expr, expr> powsum(const int* first, const int* last, double a, double b)
  {
    int minn = *std::min_element(first, last);
    double factA = 0.0;
    double factB = 0.0;
    for (; first != last; ++first)
    {
      factA += std::exp(a * (*first - minn));
      factB += std::exp(b * (*first - minn));
    }
    return{ { a * minn, factA }, { b * minn, factB } };
  }


  std::pair<expr, expr> powsum(const int* first, const int* last, memorizer const& ma, memorizer const& mb, int minn)
  {
    double factA = 0.0;
    double factB = 0.0;
    for (; first != last; ++first)
    {
      factA += ma.val(*first - minn);
      factB += mb.val(*first - minn);
    }
    return{ { ma.y() * minn, factA }, { mb.y() * minn, factB } };
  }

}
