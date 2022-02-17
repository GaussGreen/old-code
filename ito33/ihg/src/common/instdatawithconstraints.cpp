/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/instdatawithconstraints.h
// Purpose:     implementation for ihg backward instdata class
// Author:      Wang
// Created:     2004/02/13
// RCS-ID:      $Id: instdatawithconstraints.cpp,v 1.13 2006/01/27 14:26:54 yann Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/array.h"

#include "ito33/pricing/constraints.h"
#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/event.h"

#include "ihg/model.h"
#include "ihg/instdatawithconstraints.h"

using ito33::ihg::InstDataWithConstraints;

InstDataWithConstraints::InstDataWithConstraints
(
  ito33::pricing::Params& params, 
  ito33::ihg::Model& model, 
  ito33::pricing::MeshManager& meshes
) : BackwardInstData(params, model, meshes)
{  
}

void InstDataWithConstraints::Alloc(size_t nNbS)
{
  BackwardInstData::Alloc(nNbS);

  // Allocate space, and initialize to zero, since these values are
  // referenced in numoutput before timestepping begins
  m_piFrozenFlags = Array<int>(nNbS);
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    m_piFrozenFlags[nIdx] = 0;
}

void InstDataWithConstraints::UpdateBeforeStep()
{
  if(m_bHasEvent)
    for(size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
      m_piFrozenFlags[nIdx] = 0;

  BackwardInstData::UpdateBeforeStep();

}

void InstDataWithConstraints::ApplyConstraintsAfterEventsOrConstraintsUpdate()
{
  if ( m_bHasEvent )
  {
    if ( m_pConstraints )
    {
      m_pConstraints->Apply(m_pdPrices.Get(), m_piFrozenFlags.Get(), m_nNbS);

      if (m_bComputeVega)
        pricing::ApplyConstraintsToGreek
                (
                 m_pdVegas.Get(), 
                 m_piFrozenFlags.Get(), 
                 m_nNbS
                );
    }
  }
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

