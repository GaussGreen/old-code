/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Local_Functional.cpp
 *
 *  \brief this class represents the smile functionals
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2007
 */

/// remove identified warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/local_functional.h"
#include "gpmodels/Heston_Fx.h"

// gpbase
#include "gpbase/vectormanip.h"
#include "gpbase\cloneutilityfunc.h"

// gpnumlib
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/numfunction.h"

// gpcalib
#include "gpcalib/densityfunctors.h"

// gpclosedform
#include "gpclosedforms/normal.h"

#include "gpbase/eventviewerfwd.h"

#define PRICGLPOINTSNB	120
#define INVSQRT2PI		0.398942280401433


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_LocalFunctional
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_LocalFunctional
////////////////////////////////////////////////////
ARM_LocalFunctional::ARM_LocalFunctional(
		int sizeGrid,
		double nbStdDev)
:	itsFunctionals(0),
	itsGrids(0),
	itsResetTimes(0),
	itsGridSize(sizeGrid),
	itsNbStdDev(nbStdDev)
{
}


////////////////////////////////////////////////////
///	Class   : ARM_LocalFunctional
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_LocalFunctional::ARM_LocalFunctional(const ARM_LocalFunctional& rhs)
:	itsResetTimes(rhs.itsResetTimes)
{
	DuplicateCloneablePtrVectorInPlace<std::vector<double>> (rhs.itsGrids, itsGrids);
	DuplicateCloneablePtrVectorInPlace<std::vector<double>> (rhs.itsFunctionals, itsFunctionals);
}

////////////////////////////////////////////////////
///	Class   : ARM_LocalFunctional
///	Routines: toString
///	Returns : 
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_LocalFunctional::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
	size_t size = itsGrids[0]->size();
    os << indent << "Rate Functionals\n";
	size_t n=20;
	
	if (size>0)
	{
		for (size_t i=0;i<itsResetTimes.size();i++)
		{
			os << indent << "ResetTime = " << itsResetTimes[i];
			if ((i < itsFwds.size()) && (i < itsVols.size()))
				os << " Fwd = " << itsFwds[i] << " Vol = " << itsVols[i];
			os << "\n";
			os << indent << "ResetTime = " << itsResetTimes[i] << "\n";
			for (size_t k=0;k<n;k++)
			{
				os << indent << "\t" << k*(size-1)/(n-1) << "\t" << CC_NS(std,fixed) << (*itsGrids[i])[k*(size-1)/(n-1)];
				os << indent << "\t\t" << (*itsFunctionals[i])[k*(size-1)/(n-1)] << "\n";
			}
		}
	}
	os <<"\n";
	return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_LocalFunctional
///	Routines: toString
///	Returns : 
///	Action  : Object dump
////////////////////////////////////////////////////
void ARM_LocalFunctional::Calibrate(
		double resetTime,
		double fwd,
		double vol,
		ARM_DensityFunctor& density,
		bool rescaling)
{

	double mat = resetTime/K_YEAR_LEN;

	ARM_GP_VectorPtr proba = ARM_GP_VectorPtr(new ARM_GP_Vector(GetGridSize()));
	ARM_GP_VectorPtr grid = ARM_GP_VectorPtr(new ARM_GP_Vector(GetGridSize()));

	for (int k=0;k<itsGridSize;k++)
	{
		(*grid)[k]=(k-itsGridSize/2)*2.*itsNbStdDev/(itsGridSize-1.);
	}

	itsGrids.push_back(grid);

	std::vector<double>::iterator iter = proba->begin();
	std::vector<double>::iterator iter2 = grid->begin();
	for (;iter!=proba->end();++iter,++iter2)
		(*iter) = NormalCDF((*iter2));

	ARM_GP_VectorPtr func = density.Quantile(proba,fwd,mat);

	setVol(vol);
	setFwd(fwd);
	setFunc(func);
	setTime(resetTime);

	if (rescaling)
	{
		double numStateMin = grid->Elt(0);
		double numStateMax = grid->Elt(grid->size()-1);

		GaussLegendre_Coefficients glc( PRICGLPOINTSNB, numStateMin, numStateMax);
		std::vector<double> states (PRICGLPOINTSNB);
		for (size_t i(0); i<PRICGLPOINTSNB; i++)
			states[i] = glc.get_point(i);

		ARM_GP_VectorPtr rate(new std::vector<double>(PRICGLPOINTSNB));
		for (i=0; i<PRICGLPOINTSNB; i++)
			(*rate)[i] = fwd+vol*states[i];

		double state,result=0.,resultByState;
		ARM_GP_VectorPtr index = Func(resetTime,rate);
		for( i=0 ; i<PRICGLPOINTSNB; ++i )
		{
			state = glc.get_point(i);
			resultByState = index->Elt(i);
			result += glc.get_weight(i) * resultByState * exp(-0.5*state*state);
		}
		result *= INVSQRT2PI;
		double coeff = fwd/result;

		int idx = itsFunctionals.size()-1;
		for (iter=itsFunctionals[idx]->begin();iter!=itsFunctionals[idx]->end();++iter)
			(*iter) *= coeff;
	}
}

const double MINSTRIKE = 1;
const double MINPROBA = 1e-4;
const double MAXPROBA = 0.9999;

////////////////////////////////////////////////////
///	Class   : ARM_LocalFunctional
///	Routines: toString
///	Returns : 
///	Action  : Object dump
////////////////////////////////////////////////////
void ARM_LocalFunctional::CalibrateHestonFx (
		double resetTime,
		double fwd,
		const ARM_HestonModel_Fx* hestonModelFx,
		ARM_DensityFunctor& density)
{

	double mat = resetTime/K_YEAR_LEN;

	ARM_GP_VectorPtr proba = ARM_GP_VectorPtr(new ARM_GP_Vector(GetGridSize()));
	ARM_GP_VectorPtr grid = ARM_GP_VectorPtr(new ARM_GP_Vector(GetGridSize()));

	double minProb = NormalCDF(-itsNbStdDev);
	double maxProb = NormalCDF(+itsNbStdDev);

	double minHestonProba = hestonModelFx->Proba(resetTime,MINSTRIKE);
	double minVal;
	if (minHestonProba>minProb)
	{
		minVal = hestonModelFx->Quantile(resetTime,minProb);
		minVal = minVal < MINSTRIKE ? MINSTRIKE : minVal;
	}
	else
		minVal = MINSTRIKE;
	double maxVal = hestonModelFx->Quantile(resetTime,maxProb);

	for (int k=0;k<itsGridSize;k++)
	{
		(*grid)[k]=minVal+(maxVal-minVal)*k/(itsGridSize-1.0);
	}

	itsGrids.push_back(grid);

	hestonModelFx->Probas(resetTime,*grid,*proba);

	// To prevent probability greater than 1.0
	for (k = 0; k < proba->size(); ++k)
	{
		if ((*proba)[k] > maxProb)
			(*proba)[k] = maxProb;

		if ((*proba)[k] < minProb)
			(*proba)[k] = minProb;
	}

	ARM_GP_VectorPtr func = density.Quantile(proba,fwd,mat);

	setFunc(func);
	setTime(resetTime);
}


////////////////////////////////////////////////////
///	Class   : ARM_LocalFunctional
///	Routine : Func
///	Returns : ARM_VectorPtr
///	Action  : evaluation of functional
////////////////////////////////////////////////////
ARM_VectorPtr ARM_LocalFunctional::Func(double evalTime,const ARM_GP_VectorPtr& values) const
{
	if (ExistsInVector(ResetTimes(),evalTime))
	{
		size_t index = IdxFromValue(ResetTimes(),evalTime);
		size_t newSize = values->size();
		size_t gridSize = GetGridSize();
		double curNumMethState, prevNumMethState, nextNumMethState;
		double storageXmin = GetGrid(index,0);
		double storageXmax = GetGrid(index,gridSize-1);
		double storageDx = (storageXmax - storageXmin) / (gridSize-1);

		std::vector<double>& result = new std::vector<double>( newSize );

		size_t j;
		for (size_t i = 0; i<newSize; i++)
		{	
			curNumMethState = (*values)[i];
		
			/// flat extrapol if out of range
			if (curNumMethState<=storageXmin)
				result->Elt(i) =  FuncValue(index,0);
			else if (curNumMethState>=storageXmax)
				result->Elt(i) =  FuncValue(index,gridSize-1);
			else
			{
				//standard interpol
				j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
				nextNumMethState = GetGrid(index,j);
				prevNumMethState = GetGrid(index,j-1);
				
				result->Elt(i) = (  (curNumMethState - prevNumMethState) * FuncValue(index,j)
								  + (nextNumMethState - curNumMethState) * FuncValue(index,j-1) ) / (nextNumMethState - prevNumMethState);
			
			}
		}
// FIXMEFRED: mig.vc8 (30/05/2007 16:15:51):cast
		return static_cast<ARM_VectorPtr>(result);

	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model: eval not in schedule" );
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

