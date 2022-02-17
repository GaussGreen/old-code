/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 08/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file fxconvertdata.h
 *
 *  \brief
 *
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date August 2006
 */
 
#ifndef _GP_CF_FXCONVERTDATA_H
#define _GP_CF_FXCONVERTDATA_H

#include "gpbase/removenagwarning.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

#include "optimization1.h"

#include "nag.h"
#include "nage04.h"

class ARM_StdPortfolio;
class ARM_BSModel;

CC_BEGIN_NAMESPACE(ARM)

#define  DefaultTolerance    1.0e+7

///////////////////////////////////////////////////////////////
/// \class ARM_
/// \brief
///  Interface class for calibration of Fx Mkt Data
///////////////////////////////////////////////////////////////
class ARM_FXMktDataToTotemFormat
{
private:
	    ARM_StdPortfolio*  itsPortfolio;
		double itsATmVol;		
		ARM_IntVector  itsIndexDeltaPuts;
		ARM_IntVector  itsIndexDeltaCalls;
		ARM_BSModel* itsMktModel;

		ARM_GP_Vector itsUnkown;
		ARM_GP_Vector itsLowBound;
		ARM_GP_Vector itsUpBound;

		/// Nag Parameters
		double itsMaxIter;
		long   itsAlgoType;
		double itsTolerance;

		///outPut
		double itsFctObj;
		size_t itsMaturityIndex;

public:
		inline ARM_GP_Vector GetUnKnown() const {return itsUnkown;}
		inline double GetFctObjective() const {return itsFctObj;}
		inline ARM_BSModel* GetMktModel() const {return itsMktModel;}

		ARM_FXMktDataToTotemFormat(const ARM_FXMktDataToTotemFormat& rhs);
		ASSIGN_OPERATOR(ARM_FXMktDataToTotemFormat)
		~ARM_FXMktDataToTotemFormat();

		ARM_FXMktDataToTotemFormat(const ARM_StdPortfolio& portfolio,
			const ARM_GP_Vector& deltaCalls,
			const ARM_GP_Vector& deltaPuts,
			double itsATmVol,
			const ARM_BSModel& model,
			const ARM_GP_Vector& InitPoint = ARM_GP_Vector(),
			const ARM_GP_Vector& LowerBound = ARM_GP_Vector(),
			const ARM_GP_Vector& UpperBound = ARM_GP_Vector(),
			double maxIter = 50,
			double algoType = Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ,
			double Tolerance = DefaultTolerance);

		void UpdateParams(const ARM_GP_Vector& unknown);
		void WeightedSquareFunc(const ARM_GP_Vector& unknown, ARM_GP_Vector& result);
		void Derivatives(const ARM_GP_Vector& unknown,ARM_GP_Matrix& fjac);
		 
		void Calibrate();
};

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

