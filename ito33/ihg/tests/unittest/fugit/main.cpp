#include "ito33/link.h"

#include "ihg/tests/run_fugit_unit_test.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCB);

int main()
{
  return RunFugitUnitTest();
}
