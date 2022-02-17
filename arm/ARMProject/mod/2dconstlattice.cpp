/*
 * $Log
 */
#include "2dconstlattice.h"

ARM_2DConstLattice::ARM_2DConstLattice(ARM_Date startree, ARM_Date endtree,
				int nbStep, double varl, double var2)
:ARM_2DLattice(startree,  endtree,  nbStep)

{
	itsStepLen = ((endtree-startree)/K_YEAR_LEN)/nbStep;
	
}
ARM_2DConstLattice::ARM_2DConstLattice(ARM_Date startree, ARM_Date endtree,
				int nbStep)
:ARM_2DLattice(startree,  endtree,  nbStep)

{
	itsStepLen = ((endtree-startree)/K_YEAR_LEN)/nbStep;
}
void ARM_2DConstLattice::CompleteDeltaX( double MaxVar1, double MaxVar2)
{
	itsDeltaX1 = sqrt(MaxVar1);
	itsDeltaX2 = sqrt(MaxVar2);
}


