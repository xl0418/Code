/*! \file jc.h
*  \brief Entry point Janzen-Connell model
*  \author Hanno Hildenbrandt
*/

#ifndef JC_H_INCLUDED
#define JC_H_INCLUDED

#include <memory>


namespace jc {


  struct Parameter
  {
    int N = 100000;                     ///< # individuals
    double v = 0.0001;                  ///< Speciation rate
    double phi = 0.1;                   ///< Strength of phylogenetic effect
    double psi = 0.01;                  ///< Strength of abundance effect
    int64_t ticks = 10000000;           ///< number of turnovers
    std::string filename = "res.m";     ///< file name result file
    int64_t log_interval = 10;          ///< log interval for M, D and R [events]
    bool implicit = false;              ///< implicit flag
    bool verbose = false;
    int rep = 0;                        ///< repetition
  };


  std::unique_ptr<class ModelBase> CreateModel(Parameter const& param);


} // namespace jc


#endif
