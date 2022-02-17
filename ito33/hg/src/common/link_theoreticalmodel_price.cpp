/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/link_theoreticalmodel_price.cpp
// Purpose:     Force linking of theoreticalmodel_price*.cpp for HG model
// Created:     2005/01/17
// RCS-ID:      $Id: link_theoreticalmodel_price.cpp,v 1.10 2006/08/14 09:21:37 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// this file can be included in the main project when linking of all 
// theoreticalmodel_price*.cpp is required.

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceOption);

ITO33_FORCE_LINK_MODULE(HGPriceOneTouch);

ITO33_FORCE_LINK_MODULE(HGPriceCDS);

ITO33_FORCE_LINK_MODULE(HGPriceEDS);

ITO33_FORCE_LINK_MODULE(HGPriceVarianceSwap);

ITO33_FORCE_LINK_MODULE(HGPriceVarianceSwaption);

ITO33_FORCE_LINK_MODULE(HGPriceParBond);

ITO33_FORCE_LINK_MODULE(HGPriceBond);

ITO33_FORCE_LINK_MODULE(HGPriceCB);

ITO33_FORCE_LINK_MODULE(HGPriceCBOption);

ITO33_FORCE_LINK_MODULE(HGPriceReset);

ITO33_FORCE_LINK_MODULE(HGPricePEPSLike);

ITO33_FORCE_LINK_MODULE(HGPriceGeneralizedPEPSLike);

ITO33_FORCE_LINK_MODULE(HGPricePERCSLike);
