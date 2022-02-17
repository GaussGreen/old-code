#ifndef _ARMMINIMIZE_H
#define _ARMMINIMIZE_H

#include "calibration.h"

class  ARM_FRMModelMixture;

class ARM_Minimize : public ARM_Calibration
{
public:
	ARM_Minimize();

	ARM_Minimize(ARM_FRMModelMixture* model,
		const ARM_Vector& lowerBound,
		const ARM_Vector& upperBound);

	virtual ~ARM_Minimize(void);

	void Init();
	void BitwiseCopy(const ARM_Object* calib);	
	void Copy(const ARM_Object* calib);
	ARM_Object* Clone(void);
		
public:
	void StochasticMinimisation(ARM_Vector& x, int index);
	void Calibrate(ARM_Vector& x, int index);
	void sortVolParam(ARM_Vector& x);
	bool testValidSolution(const ARM_Vector& x);
private :
	ARM_Vector itsLowerBound;
	ARM_Vector itsUpperBound;
	ARM_FRMModelMixture* itsModel;
};

#endif
