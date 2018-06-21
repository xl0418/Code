/*! \file iomatlab.h
*  \brief input-output from Matlab m-files
*  \author Hanno Hildenbrandt
*/

#ifndef JANZECONNELL_IOMATLAB_H_INCLUDED
#define JANZECONNELL_IOMATLAB_H_INCLUDED

#include <fstream>
#include <utility>
#include "rndutils.hpp"
#include "jc.h"
#include "square_buffer.h"
#include "gdm.h"


namespace jc {


  // Event record
  struct event
  {
    int64_t T;                // time
    std::size_t NS;           // number of species after the event
    int sp;                   // affected species
    int ancestor;             // ancestor species
  };


  
  std::ostream& operator<<(std::ostream& os, event const& e);


  class MatLogger
  {
  public:
    MatLogger(Parameter const& param);
    ~MatLogger();

    void logEvents(std::vector<event> const& events);
    void logSnapshot(int64_t T, GDM const& D, std::vector<int> const& R);
    void logLastSnapshot(int64_t T, GDM const& D, std::vector<int> const& R);
    void logState(Parameter const& param);

  private:
    std::ofstream os_;
  };


}

#endif
