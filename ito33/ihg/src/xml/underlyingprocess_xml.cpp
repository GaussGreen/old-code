  /////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/xml/underlyingprocess_xml.cpp
// Purpose:     Restore UnderlyingProcess object from IHG XML document
// Created:     2006/06/02
// RCS-ID:      $Id: underlyingprocess_xml.cpp,v 1.2 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004-2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/vector.h"

#include "ito33/ihg/underlyingprocess.h"

#include "ihg/xml/underlyingprocess.h"
#include "ihg/xml/volatility.h"
#include "ihg/xml/hazardrate.h"

#include "ito33/xml/finance/underlyingprocess.h"
#include "ito33/xml/read.h"
#include "ito33/xml/read_vector.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

  using namespace XML;

namespace ihg
{
  

shared_ptr<UnderlyingProcess> 
XML::ReadUnderlyingProcess(const xml::node& node)
{
  shared_ptr<Volatility> pVol = ReadVolatility(node);

  shared_ptr<HazardRate> pHR = ReadHazardRate(node);

  double dPostDefaultVolatility = GetPostDefaultVolatilityFromNode(node);

  ASSERT_MSG( pVol && pHR, "Vol or hazard rate missing from the inhomogeneous "
    "underlying process" );

  shared_ptr<UnderlyingProcess>
    pUnderlyingProcess( new UnderlyingProcess(pVol, pHR) );
  
  pUnderlyingProcess->SetPostDefaultVolatility( dPostDefaultVolatility );

  return pUnderlyingProcess;
}


} // namespace hg

} // namespace ito33
