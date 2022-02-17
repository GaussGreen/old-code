/////////////////////////////////////////////////////////////////////////////
// Name:        hg/priceforwardoption.h
// Purpose:     Use HG model to price forward option
// Created:     2005/05/22
// RCS-ID:      $Id: priceforwardoption.h,v 1.3 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _HG_PRICEFORWARDOPTION_H_
#define _HG_PRICEFORWARDOPTION_H_

namespace ito33 
{ 
  namespace finance
  {
    class ForwardOption;
  }
  
namespace hg 
{

class TheoreticalModel;
class MultiOutput;

shared_ptr<MultiOutput>
PriceForwardOption(const TheoreticalModel& model,
                   const finance::ForwardOption& forwardOption);

}

}

#endif // #ifndef _HG_PRICEFORWARDOPTION_H_
