/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/transitionprobabilityoutput.h
// Purpose:     hg output class for transition probability
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobabilityoutput.h,v 1.3 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/transitionprobabilityoutput.h

    hg output class for transition probability
 */

#ifndef _ITO33_HG_TRANSITIONPROBABILITYOUTPUT_H_
#define _ITO33_HG_TRANSITIONPROBABILITYOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/hg/dlldecl.h"

namespace ito33
{

namespace finance { class ITO33_DLLDECL ModelOutput; }

namespace hg
{

/// hg output class for transition probability
class TransitionProbabilityOutput
{
public:

  void SetSpaceMesh(const std::vector<double>& pReturns)
  { 
    m_pReturns = pReturns;
  }

  const std::vector<double>& GetSpaceMesh() const
  {
    return m_pReturns;
  }

  void SetProbabilities
       (const std::vector< shared_ptr<finance::ModelOutput> >& ppProbas)
  {
    m_ppProbas = ppProbas;
  }

  const std::vector< shared_ptr<finance::ModelOutput> >&
  GetProbabilities() const 
  {
    return m_ppProbas;
  }

  void SetProbabilitiesToDefault(const std::vector<double>& pdProbaToDefaults)
  {
    m_pdProbaToDefaults = pdProbaToDefaults;
  }

  const std::vector<double>& GetProbabilitiesToDefault() const
  {
    return m_pdProbaToDefaults;
  }

private:

  std::vector<double> m_pReturns;

  std::vector< shared_ptr<finance::ModelOutput> > m_ppProbas;

  std::vector<double> m_pdProbaToDefaults;
};

} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_TRANSITIONPROBABILITYOUTPUT_H_
