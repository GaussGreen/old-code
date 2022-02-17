/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/termstructurecds_xml.cpp
// Purpose:     Restore term structure cds object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-24
// RCS-ID:      $Id: termstructure_xml.cpp,v 1.4 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructurederivative.h"
#include "ito33/finance/termstructure_enumerator.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/termstructure.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;
                           

namespace ito33
{

namespace XML
{


void ReadTermStructures
              (
                const xml::node& node,
                finance::TermStructureEnumerator& termstructures
              )
  {
    shared_ptr<finance::TermStructureCDS> ptsCDS;
    if(Restore(node, ptsCDS))
      termstructures.SetTermStructureCDS(ptsCDS);

    shared_ptr<finance::TermStructureParBond> ptsParBond;
    if(Restore(node, ptsParBond))
      termstructures.SetTermStructureParBond(ptsParBond);

    shared_ptr<finance::TermStructureEDS> ptsEDS;
    if(Restore(node, ptsEDS))
      termstructures.SetTermStructureEDS(ptsEDS);

    shared_ptr<finance::TermStructureOption> ptsOption;
    if(Restore(node, ptsOption))
      termstructures.SetTermStructureOption(ptsOption);
  }

}

}
