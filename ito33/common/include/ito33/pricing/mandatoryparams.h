/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/mandatoryparams.h
// Purpose:     Mandatory params class
// Created:     2004/08/20
// RCS-ID:      $Id: mandatoryparams.h,v 1.16 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/mandatoryparams.h
    @brief mandatory params class

    Implementation of the params class for Mandatory.
 */

#ifndef _ITO33_PRICING_MANDATORYPARAMS_H_
#define _ITO33_PRICING_MANDATORYPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/date.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/mandatory.h"

namespace ito33 
{

namespace pricing
{


/// Pricing parameters for Mandatory
class MandatoryParams : public CBLikeParams
{
public: 

  /**
      Ctor by Mandatory contract objet. 

      Other members must be initialized by the SetXXX() functions.

      @param mandatory reference to Mandatory contract objet
   */
  MandatoryParams(Mandatory& mandatory) 
                : CBLikeParams(mandatory), m_mandatory(mandatory)
  {  
  }

  /**
      Ctor by CB contract objet. This Ctor is essentailly the same
      as above, except that the memory for the contract class will
      be managed internally through the autoptr.

      Other members must be initialized by the SetXXX() functions.

      @param mandatory autoptr to Mandatory contract objet
   */
  MandatoryParams( AutoPtr<Mandatory> mandatory) 
    : CBLikeParams(*mandatory), m_mandatory(*mandatory),  
      m_clonedMandatory(mandatory)
  {  
  }

  /**
      Ctor by Mandatory contract and common financial and numerical datas. 

      @param mandatory reference to Mandatory contract objet
      @param sessionData reference to financial SessionData
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  MandatoryParams(Mandatory& mandatory,
                  const finance::SessionData& sessionData,
                  const shared_ptr<numeric::NumParams>& pNumParams,
                  const shared_ptr<numeric::MeshParams>& pMeshParams)
    : CBLikeParams(mandatory, sessionData, pNumParams, pMeshParams),
      m_mandatory(mandatory)
  {
  }

  // default dtor is ok

  /// Clone
   MandatoryParams* Clone() const;

  /// Gets a reference to the CB contract
  Mandatory& GetMandatory() const 
  { 
    return m_mandatory; 
  }
  
  /// Indicate if problem is path dependent
  virtual bool IsPathDependent() const;
  
  /// Construct path dependent grid
  virtual std::vector<double> 
  ConstructPathDepGrid(Model& model) const;

  /// Construct each path
  virtual std::vector< AutoPtr<CBLikeParams> >
  ConstructPaths(const std::vector<double>& grid) const;

  /// Construct event list
  virtual std::list< shared_ptr<PathDepEvent> >
  ConstructPathDepEvents
  (const std::vector< AutoPtr<CBLikeParams> >& paths) const;

  /// Indicate which path to save
  virtual size_t GetPathToSave(const std::vector<double>& grid) const;

  /// Initialize each path
  virtual void InitPaths(PathDepStructure& pathDepStruct);


protected:

  /// mandatory
  Mandatory& m_mandatory;

  /// Cloned object needs to manage memory of contract
  AutoPtr<Mandatory> m_clonedMandatory;

  ///indicate if averaging is supported
  bool m_bHasAveraging;

private:

  NO_COPY_CLASS(MandatoryParams);

}; // class MandatoryParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_MANDATORYPARAMS_H_
