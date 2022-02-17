/////////////////////////////////////////////////////////////////////////////
// Name:        hg/priceonetouches.h
// Purpose:     Use HG model to price multiple one touches by using homogeneity
// Created:     2006/02/24
// RCS-ID:      $Id: priceonetouches.h,v 1.2 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _HG_PRICEONETOUCHES_H_
#define _HG_PRICEONETOUCHES_H_

namespace ito33 
{ 
  namespace finance
  {
    class OneTouches;
  }
  
namespace hg 
{

class TheoreticalModel;
class MultiOutput;

shared_ptr<MultiOutput>
PriceOneTouches(const TheoreticalModel& model,
                const finance::OneTouches& oneTouches);

} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_PRICEONETOUCHES_H_
