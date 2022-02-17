/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: OUProcess.h,v $
 * Revision 1.1  2004/07/30 15:58:48  jmprie
 * Initial revision
 *
 *
 */



/*! \file OUProcess.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_OUPROCESS_H
#define _INGPMODELS_OUPROCESS_H

#include "gpbase/port.h"
#include "gpbase/curvetypedef.h"


CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \struct OUProcess
// \brief Class for functions related to Orstein-Ulenbeck process :
//        dX(t) = -mrs(t).X(t) + sigma(t)dW(t)
//-----------------------------------------------------------------------------
struct ARM_OUProcess 
{

    static double Drift(double a, double b, double mrs);

    static double Variance(double a, double b, const ARM_Curve& sigma, double mrs);

    static double IntegrateScaledSigma(double a, double b, const ARM_Curve& sigma, double scale);

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

