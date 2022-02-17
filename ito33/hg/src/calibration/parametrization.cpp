/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/calibration/parametrization.cpp
// Purpose:     implementation for parametrization class
// Created:     2005/05/20
// RCS-ID:      $Id: parametrization.cpp,v 1.25 2006/08/22 20:59:17 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/debugparameters.h"

#include "ito33/finance/error.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/eds.h"

#include "ito33/finance/basket_goodtype.h"

#include "ito33/numeric/exception.h"

#include "hg/translator.h"
#include "hg/calibratorgeneral.h"

#include "ito33/hg/underlyingprocess.h"
#include "ito33/hg/parametrization.h"
#include "ito33/hg/version.h"
#include "ito33/hg/error.h"
#include "ito33/hg/theoreticalmodel.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"

#include "hg/xml/common.h"
#include "hg/xml/parametrization.h"

extern const ito33::finance::Error 
  ITO33_CALIBRATION_FAIL,
  ITO33_EMPTY_DERIVATIVELIST,
  ITO33_INVALID_UNDERLYINGPROCESS;

extern const ito33::hg::Error ITO33_HG_CALIBRATIONFLAGS;

namespace ito33
{

namespace hg
{


Parametrization::Parametrization(const shared_ptr<UnderlyingProcess>& pUP)
                               : m_pUnderlyingProcess(pUP)                                
{
  CHECK_COND(pUP, ITO33_INVALID_UNDERLYINGPROCESS);
}

std::string Parametrization::GetDebugOutputFile() const
{
  return m_pDebug->GetDebugOutputFile("hg_calibration.xml");
}

void Parametrization::SetCalibrationFlags
                      (const std::vector<bool>& calibrationFlags)
{
  CHECK_COND( !calibrationFlags.empty(), ITO33_HG_CALIBRATIONFLAGS );

  m_calibrationFlags = calibrationFlags;
}

shared_ptr<UnderlyingProcess>
Parametrization::Calibrate(const finance::Derivatives& derivatives)
{  
  
  // Check flags and create the translator
  size_t nNbParams = m_pUnderlyingProcess->GetNbParameters();

  if ( !m_calibrationFlags.empty() ) // individual flags specified
  {
    CHECK_COND(m_calibrationFlags.size() >= nNbParams,
               ITO33_HG_CALIBRATIONFLAGS);

    m_calibrationFlags.resize(nNbParams);
  }
  else
    m_calibrationFlags.resize(nNbParams, true);
   
  Translator translator(*m_pUnderlyingProcess);
  translator.SetFlags(m_calibrationFlags);

  CHECK_COND( !derivatives.GetAll().empty() , ITO33_EMPTY_DERIVATIVELIST);

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

        const finance::Derivatives::Elements& elements( derivatives.GetAll() );
        finance::Derivatives::Elements::const_iterator iter;

        for (iter = elements.begin(); iter != elements.end(); ++iter)
        {
          tagWeights.Element(XML_TAG_DERIVATIVEWEIGHT)(iter->second);
        }
      } // weights

    } // the calibration section
  }

  CalibratorGeneral calibrator(*m_pProcess);
  
  try
  {
    m_pCalibratedUP = calibrator.Calibrate(&translator, derivatives);
  }
  catch(const ito33::numeric::Exception&)
  {
    m_pCalibratedUP = calibrator.GetLastCalibratedProcess();

    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  return m_pCalibratedUP;
}

void Parametrization::Dump(ito33::XML::Tag& tagParent) const
{
  if (m_pUnderlyingProcess)
  {
    ito33::XML::Tag 
      tagInitialGuess(XML_TAG_PARAMETRIZATION_INIT, tagParent); 

    m_pUnderlyingProcess->Dump(tagInitialGuess);
  }
}

void Parametrization::Calibrate(const finance::BasketGoodType& basket)
{
  ASSERT( basket.GetDerivatives() );

  Calibrate( *basket.GetDerivatives() );
}

shared_ptr<finance::TheoreticalModel> Parametrization::GetTheoreticalModel()
{
  shared_ptr<TheoreticalModel> 
    pTM( new TheoreticalModel(m_pCalibratedUP) );

  return pTM;
}

} // namespace hg

} // namespace ito33
