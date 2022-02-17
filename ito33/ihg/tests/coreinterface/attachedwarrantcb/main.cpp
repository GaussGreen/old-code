

#include "ito33/link.h"

#include "tests.h"

using namespace ito33;

ITO33_FORCE_LINK_MODULE(IHGPriceCB);
ITO33_FORCE_LINK_MODULE(IHGPriceAttachedWarrantConvertibleBond);

int main()
{

  TestBasicPricing();
/*
  TestShareFactor();

  TestStrike();

  TestCap();

  TestCallPutPricing();
*/

 
 return 0;
}
