/*! \file jc.h
*  \brief Entry point Janzen-Connell model
*  \author Hanno Hildenbrandt
*/

#ifndef JC_H_INCLUDED
#define JC_H_INCLUDED

#define INTERACTIVE


namespace jc {


  struct Parameter
  {
    int L = 333;                        ///< Area size
    double v = 0.00001;                 ///< Speciation rate
    double phi = 0.1;                   ///< Strength of phylogenetic effect
    double psi = 0.01;                  ///< Strength of abundance effect
    double sigmaA = 0.5;                ///< sigma abundance Gaussian
    double sigmaB = 5.0;                ///< sigma distribution Gaussian
    int64_t ticks = 10000000;           ///< number of turnovers
    std::string filename = "res.m";     ///< file name result file
    int64_t log_interval = 10;          ///< log interval for M, D and R [events]
    int64_t log_first = 0;              ///< mandatory logging of \p log first events
    double elapsed_time = 0.0;
    bool verbose = false;               ///< verbose cmd line output flag
    bool continuation = false;          ///< continuation flag
    bool dominant = false;              ///< dominant flag
    bool implicit = false;              ///< implicit flag
    bool profile = false;
    mutable int64_t totalTicks = 0;     ///< total ticks
  };


  void RunModel(Parameter const& param);


} // namespace jc


#endif
