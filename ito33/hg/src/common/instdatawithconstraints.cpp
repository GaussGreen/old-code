/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/instdatawithconstraints.cpp
// Purpose:     implementation for HG instdata for non linear system
// Created:     2004/02/13
// RCS-ID:      $Id: instdatawithconstraints.cpp,v 1.9 2006/01/24 15:13:20 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/array.h"

#include "ito33/pricing/constraints.h"
#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/event.h"

#include "hg/model.h"
#include "hg/instdatawithconstraints.h"

namespace ito33
{

namespace hg
{


InstDataWithConstraints::InstDataWithConstraints
    (pricing::Params& params, Model& model, pricing::MeshManager& meshes) 
   : BackwardInstData(params, model, meshes)
{  
}

void InstDataWithConstraints::Alloc(size_t nNbX)
{
  BackwardInstData::Alloc(nNbX);

  // The flags associated to constraints
  m_piFrozenFlags = Array<int>(nNbX);
}

void InstDataWithConstraints::UpdateBeforeStep()
{
  // Reset the flags to zero in case of an event, seems faster than use the 
  // previous which might be too wrong to be used as initial guess
  if (m_bHasEvent)
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      m_piFrozenFlags[nIdx] = 0;

  BackwardInstData::UpdateBeforeStep();
}

void InstDataWithConstraints::ApplyConstraintsAfterEventsOrConstraintsUpdate()
{
  if ( m_bHasEvent && m_pConstraints )
  {
    m_pConstraints->Apply(m_pdPrices.Get(), m_piFrozenFlags.Get(), m_nNbX);

    if ( m_sensitivityMethod == SensitivityMethod_PDE )
      for (size_t nIdxD = 0; nIdxD < m_pSensitivityData.size(); nIdxD++)
      {
        pricing::ApplyConstraintsToGreek
            (m_ppdSensitivities[nIdxD].Get(), m_piFrozenFlags.Get(), m_nNbX);
      } // loop over sensitivities

  } // if event and constraints
}

void InstDataWithConstraints::DoEvents()
{
  // No event is found and applied yet.
  m_bHasEvent = false;

  ApplyEvents();
  
  ApplyConstraintsAfterEventsOrConstraintsUpdate();
  
  const ito33::pricing::Event* pEvent;

  while ( ( pEvent = m_params.GetBasicEventAppliedAfterConstraints() ) != 0 )
    ApplyEvent(pEvent);
}


} // namespace hg

} // namespace ito33
