/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QModelAnalytics.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPMODELS_QMODELANALYTICS_H
#define _INGPMODELS_QMODELANALYTICS_H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"

/// kernel
#include <glob/armdef.h>

/// gpnumlib
#include "gpnumlib/numfunction.h"

CC_BEGIN_NAMESPACE( ARM )

/// a simple structure for the analytics of the QModel
struct QModelAnalytics
{
    private :
    static double BSFormulaForQModel(double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);

    public :
	static double BSQFunction( 
		double forward, 
		double strike, 
		double qSigma,
		double maturity,
		double qParameter,
		double df,
		int callOrPut,
        double qForward,
        double* vega=NULL);

	static double DigitalBSQFunction( 
		double forward, 
		double strike, 
		double qSigma,
		double maturity,
		double qParameter,
		double df,
		int callOrPut,
        double qForward);

    static double NormalVolToQVol(
		double forward, 
        double norSigma,
        double maturity,
        double qParameter);
};

struct QPriceFunct : public CC_NS( ARM_GP, UnaryFunc)<double,double>
{
private:
    double  itsForward;
    double  itsStrike;
    double  itsMaturity;
    double  itsQParameter;
    int     itsCallPut;
    double  itsQForward;


public:
    QPriceFunct(double F, double K, double T, double Q, double QF, int CP = K_CALL)
        : itsForward(F), itsStrike(K), itsMaturity(T), itsQParameter(Q), itsCallPut(CP), itsQForward(QF) {}

    virtual double operator () (double x) const
    { return QModelAnalytics::BSQFunction(itsForward,itsStrike,x,itsMaturity,itsQParameter,1.0,itsCallPut,itsQForward); }

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

