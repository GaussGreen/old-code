/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/readicare.h
// Purpose:     Read ICARE calibration output file
// Created:     2005/06/20
// RCS-ID:      $Id: readicare.h,v 1.2 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

namespace ito33
{

  namespace finance
  {
    class Derivatives;
  }

namespace hg
{

  class Parametrization;


void 
ReadICARE
(std::ifstream& sIn, shared_ptr<Parametrization>& pParamtrization,
 shared_ptr<finance::Derivatives>& pDerivatives);


} // namespace hg

} // namespace ito33
