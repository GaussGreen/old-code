/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization.h
// Purpose:     ihg abstract parametrization (calibration) class
// Author:      ZHANG Yunzhi
// Created:     Feb 24, 2004
// RCS-ID:      $Id: parametrization.h,v 1.30 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2003 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization.h
    @brief ihg base parametrization (calibration) class
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_H_
#define _ITO33_IHG_PARAMETRIZATION_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/parametrization.h"

#include "ito33/ihg/common.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace ihg
{

  class ParametrizationVisitor;

/**
    Base class of ihg parametrization classes.

    @noexport COM
    @nocreate
 */
class ITO33_IHG_DLLDECL Parametrization : public finance::Parametrization
{
public:

  /// ctor
  Parametrization() { }

  /// virtual dtor
  virtual ~Parametrization() { }

  virtual std::string GetDebugOutputFile() const;

  /**
      Dump all data of this parametriztion in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c for all instruments involved at once but can also be
      called "manually" if needed.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  virtual void Dump(ito33::XML::Tag& tagParent) const = 0;

  /**
      Support for visitor pattern when reading from XML.

      @param visitor the parametrization visitor object
      @noexport
   */
  virtual void Visit(ParametrizationVisitor& visitor) const = 0;

  /**
      @internal
      @brief Calibrates using a basket with the good type for instruments

      @param basket An object holding the required instruments for the 
                   calibration

      @noexport
   */
  virtual void Calibrate(const finance::BasketGoodType& basket);

  virtual shared_ptr<finance::TheoreticalModel> GetTheoreticalModel();

}; // class Parametrization


} // namespace ihg

} // namespace ito33

#endif // _ITO33_IHG_PARAMETRIZATION_H_
