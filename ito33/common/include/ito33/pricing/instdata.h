/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/instdata.h
// Purpose:     base instdata class for problem 
// Author:      Zhang Yunzhi
// Created:     2003/12/26
// RCS-ID:      $Id: instdata.h,v 1.15 2006/05/05 14:36:48 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/instdata.h
    @brief base instdata class for PDE problem depending on space variable
 */

#ifndef _ITO33_PRICING_INSTDATAT_H_
#define _ITO33_PRICING_INSTDATAT_H_

#include "ito33/array.h"

#include "ito33/pricing/instdatatimeonly.h"

namespace ito33
{

namespace pricing
{


/// base instdata class for PDE problem depending on space variable
class InstData : public InstDataTimeOnly
{
public:
  /**
     ctor
 
     @param params the Params
     @param meshes the mesh manager
   */
  InstData(Params& params, MeshManager& meshes) 
         : InstDataTimeOnly(params, meshes)
  { }

  /// dtor
  virtual ~InstData() { }


  /// @name Functions required by the Engine
 
  //@{

  // Init() : No implementation at this level, so no need to redeclare

  virtual void UpdateBeforeStep();

  // By default, do nothing at end of subgrid (reimplement in derived classes)
  virtual void UpdateAtEndOfGrid() {}
 

  // DoEvents(): Same implementation as in base, no need to reimplement

  // SetInitialValue(): No implementation at this level, so no need to redeclare

  //@}

  /// Space mesh
  const double *m_pdS;
  
  /// Log Space mesh
  const double *m_pdLogS;

  /// mesh size
  size_t m_nNbS;

  /// @name the price arrays
  //@{

  Array<double> m_pdPrices;

  Array<double> m_pdOldPrices;

  Array<double> m_pdOldOldPrices;
  
  //@}

  /// Delta array
  Array<double> m_pdDeltas;
  
  /// Gamma array
  Array<double> m_pdGammas;


protected:

  /**
     Allocate the memory for the price arrays
    
     Needs to be overloaded and called by derived class to add model specific
     items

     @param nNbS the size of the arrays to be allocated
   */
  virtual void Alloc(size_t nNbS);  

  /**
     Helper functions for swapping the price arrays
     Needs to be overloaded and called by derived class to add model specific
     items
   */
  virtual void Swap();


private:

  NO_COPY_CLASS(InstData);

}; // class InstData


} // namespace pricing

} // namespace ito33 
 
#endif // #ifndef _ITO33_PRICING_INSTDATA_H_

