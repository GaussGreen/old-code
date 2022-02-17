/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/varianceswapnumoutput.cpp
// Purpose:     implementation of NumOutput for variance swaps
// Created:     2006/03/05
// RCS-ID:      $Id: varianceswapnumoutput.cpp,v 1.4 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "hg/backwardnumoutput.h"
#include "hg/varianceswapnumoutput.h"
#include "hg/varianceswapinstdata.h"

// implement the AutoPtrDeleter for VarianceSwapNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::VarianceSwapNumOutput);
}

namespace ito33
{

  using namespace numeric;

namespace hg
{


void VarianceSwapNumOutput::Init(VarianceSwapInstData& instdata)
{  
  BackwardNumOutput::Init(instdata);
}

} // namespace hg

} // namespace ito33
