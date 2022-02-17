/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/impliedvol.cpp
// Purpose:     Class for computing implied volatility
// Author:      David
// Created:     2004/06/01
// RCS-ID:      $Id: impliedvol.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/bisecnewton.h"

#include "ihg/impliedvol.h"

namespace ito33
{

namespace ihg 
{


double ImpliedVol::Compute()
{
  double dImpliedVol = 0.0;

  shared_ptr<Volatility> pVolatility(new VolatilityFlat(0.2));

  m_theoreticalModel.SetVolatility(pVolatility);

  /*
    With a given normal hazard rate (hazard rate doesn't increase with spot), 
    the price is a strictly increasing function of the flat volatility. So this
    function will find out the flat vol if there is one in [0, 5]. Otherwise, 
    the price is either too low (< price at vol = 0) or 
    too high (> price at vol = 5)
  */

  numeric::BisecNewton solver(1.e-4);

  try
  {
    // don't begin with 5, just too high to be realistic and often problematic
    dImpliedVol = solver(*this, 0., 2);
  }
  catch(const ito33::numeric::BisecNewton::Exception& e)
  {    

    // Can have tolerance issues when the vol is close to zero. Test
    // for this by relaxing tolerance at zero boundary

    // in case that a really high price is entered that corresponding to a 
    // big vol, we check and start again  
    
    //The third argument is set to 4.99 because the max value for the volatility
    //flat is 5 then, we have to take a value strictly lower such as, even after
    //an eventual shift of the vol., this value does not become greater than 5.
    if ( fabs(solver.GetInitialFuncValue1()) < 1.e-3)
      dImpliedVol = 0.0;
    else if (solver.GetInitialFuncValue1() > 0 && solver.GetInitialFuncValue2() < 0)
      dImpliedVol = solver(*this, 2., 4.99);
    else
      throw e;
  }

  return dImpliedVol;
}


} // namespace ihg

} // namespace ito33
