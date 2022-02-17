/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/multioutput.h
// Purpose:     HG model output class for forward options
// Created:     2006/02/24
// RCS-ID:      $Id: multioutput.h,v 1.1 2006/02/24 10:06:59 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/multioutput.h
    @brief HG model output class for multiple observation points
 */

#ifndef _ITO33_HG_MULTIOUTPUT_H_
#define _ITO33_HG_MULTIOUTPUT_H_

#include "ito33/vector.h"

#include "ito33/hg/modeloutput.h"

namespace ito33
{

namespace hg
{

/// output with multiple observation points
class MultiOutput : public ModelOutput
{
public:

  MultiOutput() : ModelOutput(), m_bSensitivityOnObjectif(true) { }

  virtual ~MultiOutput() { }

  void SetPrices(const std::vector<double>& pdPrices)
  {
    m_pdPrices = pdPrices;
  }

  const std::vector<double>& GetPrices() const { return m_pdPrices; }

  void SetObjectif(double dObjectif) { m_dObjectif = dObjectif; }

  double GetObjectif() const { return m_dObjectif; }

  bool IsSensitivityOnObjectif() const { return m_bSensitivityOnObjectif; }

  void SetMultiSensitivities
  (const std::vector< std::vector<double> >& ppdSensitivities)
  {
    m_ppdSensitivities = ppdSensitivities;

    m_bSensitivityOnObjectif = false;
  }

  const std::vector< std::vector<double> >& GetMultiSensitivities() const
  {
    return m_ppdSensitivities;
  }


protected:

  std::vector<double> m_pdPrices;

  std::vector< std::vector<double> > m_ppdSensitivities;

  bool m_bSensitivityOnObjectif;

  double m_dObjectif;

private:

  NO_COPY_CLASS(MultiOutput);

};  // class MultiOutput 


} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_MULTIOUTPUT_H_
