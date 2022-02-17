#ifndef _VOLSPLINE_H
#define _VOLSPLINE_H

#include "volint.h"

class ARM_VolSplineInterpol : public ARM_VolLInterpol
{
private:
	ARM_Vector * itsDer2;
	ARM_Vector * itsVolsbyRaw;

public:


    // here m2 is not relevant
    virtual double VolatilityFunction(double m1, double k, double m2);

    virtual double VolatilityFunction(double m1, double m2);

private:

    void Init(void);

	void buildSpline();

public:

    ARM_VolSplineInterpol(void);

    ARM_VolSplineInterpol(const ARM_Date& asOf, ARM_Vector* yearTerms, 
                     ARM_Vector* strikes, ARM_Matrix* volatilities,
                     int strikeType = K_STK_TYPE_PRICE,
                     int volType = K_ATMF_VOL,
                     ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

    // Constructor for at the money vol curve 
    // (used eg for cap & floor ATM vol) 
    // -> itsStrikes contains 1 elt

    ARM_VolSplineInterpol(const ARM_Date& asOf, ARM_Vector* yearTerms, 
                     ARM_Vector* volatilities, 
                     int strikeType = K_STK_TYPE_PRICE,
                     int volType = K_ATMF_VOL,
                     ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

    ARM_VolSplineInterpol(ARM_VolFlat* volFlat, 
                     ARM_Vector* expiries = NULL,
                     ARM_Vector* undTenors = NULL);

    ARM_VolSplineInterpol(const ARM_VolSplineInterpol &);

    virtual ~ARM_VolSplineInterpol(void);

    ARM_VolSplineInterpol& operator = (const ARM_VolSplineInterpol &);

	void BitwiseCopy(const ARM_Object* srcVollint);

	virtual ARM_Object* Clone(void);

	void Copy(const ARM_Object* vollint); // !! pas hérité !!
};

#endif