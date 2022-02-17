  /////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/xml/underlyingprocess_xml.cpp
// Purpose:     Restore UnderlyingProcess object from HG XML document
// Created:     2005/04/15
// RCS-ID:      $Id: underlyingprocess_xml.cpp,v 1.3 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/vector.h"

#include "ito33/hg/underlyingprocess.h"

#include "hg/xml/underlyingprocess.h"

#include "ito33/xml/finance/underlyingprocess.h"

#include "ito33/xml/read.h"
#include "ito33/xml/read_vector.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

  using namespace XML;

namespace hg
{
  

shared_ptr<UnderlyingProcess> 
XML::ReadUnderlyingProcess(const xml::node& node)
{
  size_t nNbRegimes = GetLongFromName(node, XML_TAG_REGIME_NUMBER);

  // Resotre the volatilities
  const std::vector<double>& 
    pdVols( GetVectorFromNode<double>
            (node, XML_TAG_VOLATILITIES, XML_TAG_VOLATILITY) );

  // Restore the default intensities
  const std::vector<double>& 
    pdIntensities( GetVectorFromNode<double>
                   (node, XML_TAG_DEFAULTINTENSITIES, XML_TAG_DEFAULTINTENSITY) );

  shared_ptr<UnderlyingProcess> 
    pUnderlyingProcess( new UnderlyingProcess
                            (nNbRegimes, pdVols, pdIntensities) );

  // Restore optional jumps
  for ( xml::node::const_iterator j = node.begin(); j != node.end(); ++j )
  {
    if ( strcmp(j->get_name(), XML_TAG_JUMPS) == 0 )
    {
      size_t nIdxR1 = GetLongFromName(*j, XML_TAG_JUMP_FROM);
      size_t nIdxR2 = GetLongFromName(*j, XML_TAG_JUMP_TO);
  
      xml::node::const_iterator k;
      Jumps jumps;
      for ( k = j->begin(); k !=j->end(); ++k)
      {
        if ( strcmp(k->get_name(), XML_TAG_JUMP) == 0 )
        {
          double dIntensity = GetDoubleFromName(*k, XML_TAG_INTENSITY);
          double dAmplitude = GetDoubleFromName(*k, XML_TAG_AMPLITUDE);

          jumps.push_back(Jump(dIntensity, dAmplitude));
        }
      }

      pUnderlyingProcess->SetJumps(nIdxR1, nIdxR2, jumps);
    }
  }

  // Restore the post default volatility
  double dPostDefaultVolatility = GetPostDefaultVolatilityFromNode(node);  
  pUnderlyingProcess->SetPostDefaultVolatility( dPostDefaultVolatility );
  
  return pUnderlyingProcess;
}


} // namespace hg

} // namespace ito33
