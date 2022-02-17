/////////////////////////////////////////////////////////////////////////////
// Name:        cb_base_xml.cpp
// Purpose:     Restore Bond object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-12-05
// RCS-ID:      $Id: cb_base_xml.cpp,v 1.7 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/cb_base.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/cb_base.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/putschedule.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

void ito33::XML::GetOptionalCBBaseDataFromNode
        (const xml::node& node, CBBase& cbBase)
{
  GetOptionalConvertibleLikeDataFromNode(node, cbBase);

  shared_ptr<finance::CallSchedule> pCalls;
  if ( Restore(node, pCalls) )
    cbBase.SetCallSchedule(pCalls);

  shared_ptr<finance::PutSchedule> pPuts;
  if ( Restore(node, pPuts) )
    cbBase.SetPutSchedule(pPuts);
}
