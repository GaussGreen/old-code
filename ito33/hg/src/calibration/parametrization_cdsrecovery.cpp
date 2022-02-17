/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/calibration/parametrization_cdsrecovery.cpp
// Purpose:     implementation for parametrization class with cds recovery 
// Created:     2005/07/27
// RCS-ID:      $Id: parametrization_cdsrecovery.cpp,v 1.11 2006/08/23 15:57:05 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/sessionData.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/eds.h"

#include "ito33/numeric/exception.h"

#include "hg/translator.h"
#include "hg/calibrator_cdsrecovery.h"

#include "ito33/hg/underlyingprocess.h"
#include "ito33/hg/version.h"
#include "ito33/hg/error.h"
#include "ito33/hg/parametrization_cdsrecovery.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"

#include "hg/xml/common.h"
#include "hg/xml/parametrization.h"

extern const ito33::finance::Error ITO33_INVALID_UNDERLYINGPROCESS,
                                   ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace hg
{


ParametrizationCDSRecovery::ParametrizationCDSRecovery(shared_ptr<UnderlyingProcess> pUP)
                               : Parametrization(pUP), m_bUseSpreads(true)
{
}

shared_ptr<UnderlyingProcess>
ParametrizationCDSRecovery::Calibrate(const finance::Derivatives& derivatives)
{  

  if ( IsDebugOutputEnabled() )
  {    
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_HG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_ROOT_VERSION, ITO33_HG_VERSION_DOT_STRING);

    const finance::Derivatives::Elements& elements = derivatives.GetAll();

    derivatives.GetSessionData()->Dump(tagRoot);  

    {
      // dump the initial underlying process
      ito33::XML::Tag tagParam(XML_TAG_HG_CALIBRATION, tagRoot);      
      Dump(tagParam);
      
      // dump the derivatives to be calibrated
      {
        ito33::XML::Tag tagDerivatives(XML_TAG_DERIVATIVES, tagParam);
      
        finance::Derivatives::Elements::const_iterator iter;
      
        for (iter = elements.begin(); iter != elements.end(); ++iter)
          iter->first->Dump(tagDerivatives);
      }

      // dump the derivative weights
      {
        ito33::XML::Tag tagWeights(XML_TAG_DERIVATIVEWEIGHTS, tagParam);

        finance::Derivatives::Elements::const_iterator iter;
      
        for (iter = elements.begin(); iter != elements.end(); ++iter)
        {
          tagWeights.Element(XML_TAG_DERIVATIVEWEIGHT)(iter->second);
        }
      } // weights

    } // the calibration section
  }

  Translator translator(*m_pUnderlyingProcess);
  
  // do we support partial calibration here?

  // extra one for cds recovery value
  size_t nNbParams = m_pUnderlyingProcess->GetNbParameters() + 1;

  m_calibrationFlags.resize(nNbParams, true);
  
  CalibratorCDSRecovery calibrator(*m_pProcess);

  calibrator.CalibrateWithSpreads(m_bUseSpreads);
  
  try
  {
    m_pCalibratedUP = calibrator.Calibrate(&translator, derivatives);
    m_dCDSRecovery = calibrator.GetCDSRecovery();
  }
  catch(ito33::numeric::Exception)
  {
    m_pCalibratedUP = calibrator.GetLastCalibratedProcess();
    m_dCDSRecovery = calibrator.GetCDSRecovery();

    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  return m_pCalibratedUP;
}


} // namespace hg

} // namespace ito33
