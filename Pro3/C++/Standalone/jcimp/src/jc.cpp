#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include "rndutils.hpp"
#include <numeric>
#include <algorithm>
#include <chrono>
#include "jc_base.h"


namespace jc {


  rndutils::xorshift128 thread_local RndEng = rndutils::make_random_engine();


  extern ModelBase* CreateNeutralModel(Parameter const&);
//  extern ModelBase* CreateNeutralDModel(Parameter const&);
//  extern ModelBase* CreateNeutralRModel(Parameter const&);
  extern ModelBase* CreateImplicitModel(Parameter const&);


  std::unique_ptr<ModelBase> CreateModel(Parameter const& param)
  {
    auto p = param;
    auto pp = p.filename.find_last_of('.');
    std::string filename(param.filename.begin(), param.filename.begin() + pp);
    std::stringstream ss;
    ss << filename << '_' << param.rep << '_' << static_cast<int>(param.phi * 100000) << '_' << static_cast<int>(param.psi * 100000) << std::string(param.filename.begin() + pp, param.filename.end());
    p.filename = ss.str();
    std::unique_ptr<ModelBase> pm;
    if (p.psi == 0.0 && p.phi == 0.0) pm.reset(CreateNeutralModel(p));
    //else if (p.phi == 0.0) pm.reset(CreateNeutralDModel(p));
    //else if (p.psi == 0.0) pm.reset(CreateNeutralRModel(p));
    else pm.reset(CreateImplicitModel(p));
    return std::move(pm);
  }


} // namespace jc
