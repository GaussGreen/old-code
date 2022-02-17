//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : EntropyExpLossInterpolator.cpp
//
//   Description : Interpolates expected losses using entropy minimisation
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Potential.hpp"
#include "edginc/UniformDensityPrior.hpp"
#include "edginc/EntropyExpLossInterpolator.hpp"
#include "edginc/GridFactory.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** minimum separation between tranches */
const double EntropyExpLossInterpolator::minSep = 10*1.0e-10; 

/** distance between strikes in homogeneous part of grid */
const double EntropyExpLossInterpolator::elStep = 0.01;

/** Private constructor */
EntropyExpLossInterpolator::EntropyExpLossInterpolator() :	CObject(TYPE),
															strikes(0),
															prior(new UniformDensityPrior(0.0,1.0)),
															scale(1.0),
															optimizer1(Simplex::create(10.0)),
															optimizer2(Simplex::create(1.0)),
															showPrior(false)
{	
}

/** Public constructor */
EntropyExpLossInterpolator::EntropyExpLossInterpolator(const DoubleArray & strikes,
													   IDensityPriorSP prior) :	
															CObject(TYPE),
															strikes(new DoubleArray(strikes)),
															prior(prior),
															scale(1.0),
															optimizer1(Simplex::create(10.0)),
															optimizer2(Simplex::create(1.0)),
															showPrior(false)
{}

/** Destructor */
EntropyExpLossInterpolator::~EntropyExpLossInterpolator()
{}

////////////////////////////////////////////////////////////////////////////////
void EntropyExpLossInterpolator::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(EntropyExpLossInterpolator, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IExpectedLossInterpolator);
    EMPTY_SHELL_METHOD(defaultEntropyExpLossInterpolator);

	FIELD(strikes, "Base strikes at which to interpolate [default is 1% grid up to highest strike < 1-RR]");
	FIELD_MAKE_OPTIONAL(strikes)

	FIELD(prior, "Loss probability density prior [default: UniformDensityPrior]");
	FIELD_MAKE_OPTIONAL(prior);

	FIELD(scale ,"Postitive scale factor for expected losses [default = 1.0]");
	FIELD_MAKE_OPTIONAL(scale);

	FIELD(optimizer1, "First optimizer. [default = simplex(10)]");
	FIELD_MAKE_OPTIONAL(optimizer1);

	FIELD(optimizer2, "Second optimizer. [default = simplex(1)]");
	FIELD_MAKE_OPTIONAL(optimizer2);

	FIELD(showPrior, "TRUE if want to have loss density to be given by unmodified prior [default = FALSE]");
	FIELD_MAKE_OPTIONAL(showPrior);

}


IObject* EntropyExpLossInterpolator::defaultEntropyExpLossInterpolator() {
    return new EntropyExpLossInterpolator();
}

CClassConstSP const EntropyExpLossInterpolator::TYPE = 
    CClass::registerClassLoadMethod("EntropyExpLossInterpolator", 
                                    typeid(EntropyExpLossInterpolator), 
                                    load);

/////////////////////////////////////////////////////////////////////////////

/** Return interpolated expected loss surface obtained by
	entropy minimisation. */

//////////////////////////////////////////////////////////////////////////////
ExpectedLossSurfaceSP EntropyExpLossInterpolator::getELSurface(
		const ExpectedLossSurface & targetExpLosses	/** Expected loss points to interpolate */
		) const
{
	static const string method ="EntropyExpLossInterpolator::getELSurface";
	try
	{

		int ti, si; 
		
		/** recovery rate */
		double RR = targetExpLosses.getRecovery();

		
		/** dates at which we interpolate */
		DateTimeArrayConstSP dates = targetExpLosses.getDates();

			
		/** set range of definition for prior */
		if(prior.get() == NULL) throw ModelException("prior is NUll");

		prior->setRange(0,1-RR);

		

		/** baseStrikes for target EL's. Need to be <= 1-RR */
		DoubleArray baseStrikes;

		/** strikes corresponding to input targetExpLosses */
		DoubleArray targetStrikes = *targetExpLosses.getStrikes();
		
		for(si = 0;si < targetStrikes.size();si++)
		{
			double strike = targetStrikes[si];
			if(strike <= 1-RR) baseStrikes.push_back(strike);
		}
		
		//check that last strike is 1-RR 
		//TODO check that we need this
		//if(baseStrikes.back() != 1-RR) throw ModelException("strikes in targetExpLosses need to include 1-RR");
		
		/** number of baseStrikes */
		int nBaseStrikes = baseStrikes.size();
		if(nBaseStrikes < 1) throw ModelException("nBaseStrikes<1");

		DoubleArray dummy(0);
		
		/** strikes for final ExpectedLossSurface */
		DoubleArraySP storageStrikes = createStorageStrikes(baseStrikes,RR);

		

		/** number of strikes */
		int nStrikes = storageStrikes->size();

		/** Calculation grid: these are the grid points at which we want to evaluate the interpolated EL's
			using entropy. These are storageStrikes <= 1-RR.
			Maximum value is 1-RR since we have flat EL curve after that (which entropy can't handle). */
		DoubleArray calculationGrid(0);
		for(si = 0; si< nStrikes;si++)
		{
			if((*storageStrikes)[si] <= 1-RR) calculationGrid.push_back((*storageStrikes)[si]);
		}

		
		
		/** target base expected losses populated later */
		DoubleArray baseEL(nBaseStrikes);

		/** size of calculation grid */
		int nCalcGrid = calculationGrid.size();

		/** interpolated base losses which are populated by interpolation routine */
		DoubleArray baseELout(nCalcGrid);

		


		
		// caluclate expected loss matrix ---------------------------------------------------------------	

		// initialise expected losses 
		DoubleArrayArraySP expLosses(new DoubleArrayArray(storageStrikes->size()));
		int datesSize = dates->size();
		for(si=0;si<storageStrikes->size();si++) (*expLosses)[si]= DoubleArray(datesSize);
		
		// special case for first point in dates
		// this is probably the value date in which case all expected losses are 0

	
		// check if all expected losses at first date are 0
		if(dates->size()<1) throw ModelException("dates array is empty");
		int tStart = 0;
		int index = 0;

		/** losses in input targetExpLosses */
		DoubleArrayArray targetLosses = *targetExpLosses.getLosses();

		while(index < nBaseStrikes && targetLosses[index][tStart] < minSep) index++;
		
		// if alll expected losses are 0 then set output losses to 0
		if(index >= nBaseStrikes)
		{
			for(si = 0;si<nStrikes;si++) (*expLosses)[si][tStart] = 0.0;  
			tStart++;
		}

		
		if(scale < 0)
		{
			throw ModelException("scale "+Format::toString(scale)+" < 0");
		}
        
        // previous expected losses
        double baseELprev = 0;
        double baseELprevPrev = 0;

		// loop through all points in dates and do strike interpolation
		for(ti = tStart;ti<dates->size();ti++) 
		{
			/** prior Grid : this is grid which defines the prior between 0 and 1-RR. 
			Density is assumed piecewise const between these points */
			DoubleArraySP priorGrid = prior->getGridPoints(0, 1-RR, (*dates)[ti]);
			
			/** integration grid: These are the points at which we integrate the potential in the entropy algorithm.
			The requirement is that both the prior function and the constraint coefficent functional 
			is piecewise linear between these points. This is the merger of priorGrid and calculationGrid 
			(note that baseStrikes are already included in calculation grid  */
			DoubleArraySP integrationGrid = GridFactory::mergeNoDuplicate(calculationGrid, *priorGrid);

			// Double check that final point is less than 1-RR
			if(integrationGrid->back() > 1-RR) throw ModelException("intergrationGrid cannot have points > 1-RR"); // This can't happen

			/** coeefficient matrix needed for strike interpolation. */
			DoubleArrayArraySP coeffMatrix = constraintCoeffMatrix(
				*integrationGrid,
				baseStrikes			
				);

			// populate target expected losses
			for(si=0; si< nBaseStrikes;si++)
			{
                if(scale > 1.0)
                {
                    // baseEL is bounded by fact that we require minSep separation between scaled tranches
                    double maxEL = baseStrikes[si]; // base case
                    if(si == 1)
                    {
                        maxEL = baseELprev + (baseStrikes[si] - baseStrikes[si-1])*(baseELprev/baseStrikes[si-1] - 2*minSep);
                    }
                    if(si > 1)
                    {
                        maxEL = baseELprev;
                        maxEL += (baseStrikes[si] - baseStrikes[si-1])*((baseELprev-baseELprevPrev)/(baseStrikes[si-1]-baseStrikes[si-2]) - 2*minSep);
                    }
                
                    baseEL[si] = Maths::min(scale*targetExpLosses.getBaseEL(baseStrikes[si],(*dates)[ti]), maxEL);
                    
                    baseELprevPrev = baseELprev;
                    baseELprev = baseEL[si];


                }
                else
                {
                    baseEL[si] = scale*targetExpLosses.getBaseEL(baseStrikes[si],(*dates)[ti]);
                }
			}


			// call interpolation routine
			interpolate(
				baseStrikes,
				baseEL,
				*integrationGrid,
				*prior->getDensity(*integrationGrid,(*dates)[ti]),
				*coeffMatrix.get(),
				calculationGrid,		
				baseELout
				);
			

			
			// copy to expLosses
			for(si = 0; si<nStrikes;si++)
			{
				if(si < nCalcGrid)
				{
					(*expLosses)[si][ti] = baseELout[si];
				}
				else (*expLosses)[si][ti] = baseELout.back();
			} // end copy

		} // end loop over dates



		// copy dates and strikes to (unconst) SP's so that can create ExpectedLossSurface
		DateTimeArraySP datesSP(new DateTimeArray(*dates));
		
		return ExpectedLossSurfaceSP(new ExpectedLossSurface(datesSP,storageStrikes,expLosses,RR));
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}

////////////////////////////////////////////////////////////////////////////////////////
/** create strike at which expected losses for finale ExpectedLossSurface are defined  
	This includs 0, 1 and the baseStrikes
	*/
///////////////////////////////////////////////////////////////////////////////////
DoubleArraySP EntropyExpLossInterpolator::createStorageStrikes(
	const DoubleArray & baseStrikes, // base strikes of traget tranches
	double RR						// recovery
	) const
{
	static const string method = "EntropyExpLossInterpolator::createStrikes";
	try
	{
	
		int i;
		/** output grid */
		DoubleArraySP outGrid;

		int numStrikes = baseStrikes.size();
   
        if(numStrikes < 1) throw ModelException("Need at least one base strike");
		 
		// if strikes grid has not been populated then create default grid
		if(strikes.get() == NULL) 
		{
			outGrid = DoubleArraySP(new DoubleArray(0));
			/** last strike for homogeneous part this is last baseStrike <= 1-RR*/
			double highStrike = Maths::min(1-RR,baseStrikes[0]); // must have at least one baseStrike

			if(baseStrikes.back() < 1)
			{       
				if(numStrikes>1) highStrike = Maths::min(1-RR, baseStrikes[numStrikes-1]);
			}
			else
			{       
				if(numStrikes> 2) highStrike = Maths::min(1-RR, baseStrikes[numStrikes-2]); 
			}
			

			// start at 0
			outGrid->push_back(0.0);
			
			// populate remaining strikes
			while(outGrid->back() + elStep < highStrike)
			{
				outGrid->push_back(outGrid->back() + elStep);
			}
			
			// add highStrike is not included
			if(outGrid->back() < highStrike) outGrid->push_back(highStrike);

		}
		else // just sort the strikes array
		{	
			outGrid = DoubleArraySP(new DoubleArray(*strikes));
			sort(outGrid->begin(), outGrid->end());

			// check for duplicates
			for(i = 1;i<outGrid->size();i++) 
			{
				if((*outGrid)[i] == (*outGrid)[i-1]) throw ModelException("strikes array contains duplicate strikes ("+
																				Format::toString((*outGrid)[i])+")");
			}
					

		}
		
		// add base strikes if not included	
		outGrid = GridFactory::mergeNoDuplicate(*outGrid, baseStrikes);
		
		// add 0 , 1-RR and 1  and b if not included
		DoubleArray addStrikes(0);
		addStrikes.push_back(0.0);
		if(addStrikes.back() < 1-RR) addStrikes.push_back(1-RR);
		if(addStrikes.back() < 1.0) addStrikes.push_back(1.0);

		outGrid = GridFactory::mergeNoDuplicate(*outGrid, addStrikes);


		return outGrid;

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


/////////////////////////////////////////////////////////////////////////////////
/**
  coefficient matrix h_n(x) needed in for strike interpolation

   For the constrained minimisation we need to write the constraints in the form:
	\int_0^1 h_n(x)density(x)dx = c_n
	where h_n(x) = min(x, strike_n) and c_n = baseEL_n
*/
//////////////////////////////////////////////////////////////////////////////////////
DoubleArrayArraySP EntropyExpLossInterpolator::constraintCoeffMatrix(
	const DoubleArray & integrationGrid,
	const DoubleArray & baseStrikes
	) const
{
	static const string method = "EntropyExpLossInterpolator::entropyStrikeCoeffMatrix";
	try
	{
		int si,i;		
		
		int numStrikes = baseStrikes.size();
		int gridSize = integrationGrid.size();
		
	
		/** constraint coeef matrix */
		DoubleArrayArraySP h(new DoubleArrayArray(numStrikes));
	
		for(si = 0;si<numStrikes;si++)	
		{
			(*h)[si] = DoubleArray(gridSize);
			for(i = 0;i< gridSize;i++)
			{
				(*h)[si][i] = Maths::min(integrationGrid[i],baseStrikes[si]);  	
			}
		}
		
		return h;

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
/** 
  expected loss interpolation accross using entropy
  Inputs are array of base expected losses at base strikes and a prior distribtution
  Output is Array of base expected losses at calculation grid points
  
*/
///////////////////////////////////////////////////////////////////////////////////////////
void EntropyExpLossInterpolator::interpolate(
		const DoubleArray & baseStrikes,			// (I) strikes for target base expected losses
		const DoubleArray & baseEL,					// (I) target base expected loses. 
		const DoubleArray & integrationGrid,		// (I) integration grid. Assume h_coeffs and prior are piecewise linear between points 
		const DoubleArray & prior,					// (I) prior probability density defined at integration grid points. Assumed piecewise linear inbetweem
		const DoubleArrayArray & h_coeffs,			// (I) constraint coefficient matrix returned by constraintCoeffMatrix()
		const DoubleArray & calculationGrid,		// (I) points at which we want output losses	
		DoubleArray & baseELout						// (O) output expected losses
		) const
{
	static const string method = "EntropyTrancheQuoteTimeInterpolation::interpolate";
	try
	{
		int i,x;
		
		int integGridSize = integrationGrid.size();
		int calcGridSize = calculationGrid.size();
		int numStrikes = baseStrikes.size();
	

		// validation ----------------------------------------------------------
		if(baseELout.size() != calculationGrid.size())
		{
			throw ModelException("baseELout has incorrect strike dimension");
		}

		if(prior.size() != integGridSize)
		{
			throw ModelException("prior does not have same size as integrationGrid");
		}
		
	
		if(baseEL.size() != numStrikes) throw ModelException("baseEL has differnt size to baseStrikes");
		
		

		// check concavity. I.e. check normalised tranches have separation >= minSep
		if(baseStrikes[0] <=0) throw ModelException("First base strike needs to be > 0");
		
		double prevTranche = baseEL[0]/baseStrikes[0];
		double tranche = prevTranche;
		
		for(i = 1; i<numStrikes;i++)
		{
			if(baseStrikes[i] <= baseStrikes[i-1]) throw ModelException("baseStrikes need to be strictly increasing");

			tranche = (baseEL[i]-baseEL[i-1])/(baseStrikes[i]-baseStrikes[i-1]);
			
			if(prevTranche < tranche + minSep)
			{
				throw ModelException("EL's not concave at index: "+ Format::toString(i));
			}

			prevTranche = tranche;
		}
		// check that final tranche >= minSep
		if(tranche < minSep)
		{
			throw ModelException("final normalised tranche is less than minSep = "+ Format::toString(minSep));
		}


		// set up potential ----------------------------------------------------------------------------
		Potential phi(baseEL, h_coeffs,integrationGrid, prior, numStrikes); 
		
			
		// get lagrange multipliers that minimise the potential ------------------------------------------
		DoubleArray lagrangeMultipliers(numStrikes,0);
		
		/** first guess for lagrange multipliers: all multipliers = 0 */
		DoubleArray guess(numStrikes, 0.0);
	
		
		
		// solve using two optimizers sequentially
		// This can improve convergence especially for simplex
		// if only want to show prior then lagrange multipliers are trivial and remain 0
		if(!showPrior)
		{
			optimizer1->minimize(phi,guess,lagrangeMultipliers);
			guess = lagrangeMultipliers; 

			optimizer2->minimize(phi,guess,lagrangeMultipliers);
		}


		/*
		// start with best guess and do simplex minimisation
		SimplexSP simplex = Simplex::create(10.0);

		//guess = lagrangeMultipliers; 
		simplex->minimize(phi,guess,lagrangeMultipliers);
		
		// repeat simplex to improve convergence
		SimplexSP simplex2 = Simplex::create(1.0);
		guess = lagrangeMultipliers; 
		simplex2->minimize(phi,guess,lagrangeMultipliers);
		*/

		/** normalisation factor for density */
		double norm = phi.norm(lagrangeMultipliers);
		

 
		//-------------------------------------------------------------------------------------
		// populate losses matrix : loss given by \int_0^1 min(x,strike) \density_B(x) dx 
		// this is integral of piecewise linear function f(x) = min(x,strike) against phi which can be done exactly
		//
		// TODO: We can use EL(k+dk) = EL(k) + \int_k^1 min(max(0,x-k),dk) dx 
		
		
		/** testFunction f(x) */
		DoubleArray testFunc(integGridSize);

		/* output integral needed for integration */
		DoubleArray outIntegral(integGridSize);

		
		// loop through strorage points that we need to populate
		for(i = 0;i < calcGridSize;i++)
		{
			if(calculationGrid[i] == 0) // 0 strike case
			{
				baseELout[i] = 0.0;
			}
			else
			{
				// define function f(x) = min(x,strike) to integrate against
				for(x = 0; x< integGridSize;x++)
				{
					testFunc[x] = Maths::min(integrationGrid[x], calculationGrid[i]);
				}
			
				phi.integral_PhiXfX(
					lagrangeMultipliers,
					testFunc,
					norm,
					outIntegral
					);

				baseELout[i] = outIntegral.back();
			}
		}	


	} catch(exception& e){
		throw ModelException(e,method);
	}
}




/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool EntropyExpLossInterpolatorLoad() {
    return (EntropyExpLossInterpolator::TYPE != 0);
}


DRLIB_END_NAMESPACE


