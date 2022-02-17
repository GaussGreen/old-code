
#include "ito33/finance/computationalflags.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/pricingreader.h"

#include "fugit_test_case.h"

using namespace ito33;
using namespace ito33::finance;


void FugitTestCase::ZeroBeingExercisedCase1()
{
  ihg::XML::PricingReader reader("ihgcb_fugit_should_be_0.xml");

  shared_ptr<Derivative> pDerivative;
  reader.ReadDerivative(pDerivative);

  shared_ptr<ito33::ihg::TheoreticalModel> 
    pModel(new ito33::ihg::TheoreticalModel);
  reader.ReadTheoreticalModel(pModel);

  shared_ptr<ComputationalFlags> pFlag(new ComputationalFlags);
  pFlag->SetComputeFugit(true);
  pDerivative->SetComputationalFlags(pFlag);

  CPPUNIT_ASSERT( pModel->Compute(*pDerivative)->GetFugit() == 0 );
}

