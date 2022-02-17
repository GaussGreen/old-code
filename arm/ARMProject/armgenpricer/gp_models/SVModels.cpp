/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV.cpp
 *
 *  \brief file for the Markovian Stochastic volatility model
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */



/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/SVModels.h"
#include "gpmodels/ModelParamsMSV.h"
#include "gpmodels/riccatimsv.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/datestrip.h"

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"


#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/riccati.h"


/// gpnumlib
#include "gpnumlib/gaussiananalytics.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"
#include "gpnumlib/rungekutta.h"

/// kernel 
#include <inst/irindex.h>

CC_BEGIN_NAMESPACE( ARM )

#define FIRST_STATE_VARIABLE    0
const double VOL_LIMIT			  = 0.000001;
const double GL_NBPOINTS_FIRST    = 4; // Number of points per Year
const double GL_NBPOINTS_SECOND	  = 4; // Number of points per Year
const double ODE_PRECISION		  = 1e-5; 


////////////////////////////////////////////////////
///	Class  : ARM_GeneralSwaptionData
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_GeneralSwaptionData::ARM_GeneralSwaptionData (int sizeVector)
{
	itsBetaPart_1 = std::vector<double> (sizeVector + 1, 0.0);
	itsBetaPart_2 = std::vector<double> (sizeVector + 1, 0.0);
	itsEpsPart_1 = std::vector<double> (sizeVector + 1, 0.0);
	itsEpsPart_2 = std::vector<double> (sizeVector + 1, 0.0);
	itsEpsPart_3 = std::vector<double> (sizeVector + 1, 0.0);

}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SVModels::ARM_SVModels( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params )
:	ARM_PricingModelIR(zc,params) , 	itsIntegratedVol(0)

{
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SVModels::ARM_SVModels(const ARM_SVModels& rhs)
: ARM_PricingModelIR(rhs)
{
	itsIntegratedVol = rhs.itsIntegratedVol ;
    // Copy class attributes if any
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SVModels::~ARM_SVModels()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_SVModels& ARM_SVModels::operator=(const ARM_SVModels& rhs)
{
	itsIntegratedVol = rhs.itsIntegratedVol ;
	if(this != &rhs)
	{
		ARM_PricingModelIR::operator=(rhs);
		// Copy class attributes if any
	}
	return *this;
}

//////////////////////////////////////////////////// ////////////////////////////////////////////				
///	Class  : ARM_SVModels
///	Routine: PrecomputeDatas
///	Returns: 
///	Action : Precompute Datas to accelerate Shift , VolOfVol and Volatility computation
////////////////////////////////////////////////////////////////////////////////////////////////

double computeEpsPart (double t1, double t2, double volMrs, double volOfVol)
{
#if defined(__GP_STRICT_VALIDATION)
	if(volMrs<K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Vol Mean Reversion is Null!");
#endif
	double volMrs_2 = 2 * volMrs;
	return (volOfVol*volOfVol*(exp(volMrs_2 * t2) - exp(volMrs_2 * t1))/volMrs_2);
}

									////////////////////////////////////

void ARM_SVModels::PrecomputeDatas(
				const string& curveName,
				double resetTime, 
				double floatStartTime,
				double floatEndTime,
				const std::vector<double>& fixPayTimes,
				const std::vector<double>& fixPayPeriods)								   
{
	std::vector<double> Schedule = ((ARM_ModelParamsMSV*) GetModelParams())->GetSchedule();
	double startTime = 0.0;
	double endTime,weight_in,Volt_in,t_in,temp_in,volMrs ;
	
	
	std::vector<double> newSched;
	int size = 0;
	while((size<Schedule.size()) && (Schedule[size]< resetTime))
	{
		newSched.push_back(Schedule[size]);
		size ++;
	}
	newSched.push_back(resetTime);
	
	int sizeVector = newSched.size() + 1 ;
	itsStoredData.itsBetaPart_1 = std::vector<double> (sizeVector, 0.0);
	itsStoredData.itsBetaPart_2 = std::vector<double> (sizeVector, 0.0);
	itsStoredData.itsEpsPart_1 = std::vector<double> (sizeVector, 0.0);
	itsStoredData.itsEpsPart_2 = std::vector<double> (sizeVector, 0.0);
	itsStoredData.itsEpsPart_3 = std::vector<double> (sizeVector, 0.0);


	GaussLegendre_Coefficients c_in(GL_NBPOINTS_SECOND);
	
	for (size_t i=1; i<(newSched.size()+1) ;++i)
	{
		endTime = newSched[i-1];
		double volOfVol = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(endTime);
		endTime/=K_YEAR_LEN;

		/// We have constant Mean Reversion
		volMrs		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetValueAtPoint(0);
		itsStoredData.itsEpsPart_1[i] = itsStoredData.itsEpsPart_1[i-1] + computeEpsPart(startTime ,endTime, volMrs, volOfVol);
		itsStoredData.itsEpsPart_2[i] = itsStoredData.itsEpsPart_2[i-1] + computeEpsPart(startTime ,endTime, volMrs, 1.0);

		////////// Computation of Part 3 /////////////////////////////////

		GaussLegendre_Coefficients coeff_in(&c_in,startTime,endTime); 
		for(int m=0;m<GL_NBPOINTS_SECOND;++m)
		{
			weight_in = coeff_in.get_weight(m);
			t_in = coeff_in.get_point(m);
			Volt_in = ComputeEquivalentVolatilityt(curveName,t_in * K_YEAR_LEN,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods);
			temp_in = Volt_in*Volt_in*weight_in;			

			itsStoredData.itsBetaPart_1[i] += temp_in;
			itsStoredData.itsBetaPart_2[i] += temp_in*sinh(t_in*volMrs)/volMrs;
			itsStoredData.itsEpsPart_3 [i-1] += temp_in * exp(-volMrs * t_in);
		}
		itsStoredData.itsBetaPart_1[i]+=itsStoredData.itsBetaPart_1[i-1];
		itsStoredData.itsBetaPart_2[i]+=itsStoredData.itsBetaPart_2[i-1];
		startTime = endTime;
	}
	
	for (i=newSched.size(); i>0 ;i--)
	{
		itsStoredData.itsEpsPart_3[i-1]	 += itsStoredData.itsEpsPart_3[i];
	}
	itsStoredData.itsNewSchedule = newSched;
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: GetFloatDateStrip
///	Returns: a datestrip according to the floating convention
///	Action : computes the datestrip corresponding to the index and the 
///			default convention!
////////////////////////////////////////////////////
ARM_DateStrip* ARM_SVModels::GetFloatDateStrip( const ARM_Date& startDate, const ARM_Date& endDate,const ARM_IRIndex* pIndex ) const
{	
	//// AAAAAAAAAAAAARRRRRRRRRGGGGGGGGGGG ugly const_cast
	ARM_Currency* pCcy	= const_cast<ARM_SVModels* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatResetFreq	= pIndex->GetResetFrequency();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();
	char* floatResetCal = pCcy->GetResetCalName(pIndex->GetIndexType());
	int fwdRule			= pCcy->GetFwdRule();
	int intRule			= pIndex->GetIntRule();
	int stubRule		= K_SHORTSTART;
	int resetGap		= pIndex->GetResetGap();
	int floatPayFreq	= pIndex->GetPayFrequency();
	int floatPayGap		= pIndex->GetPayGap();
	char* floatPayCal	= pCcy->GetPayCalName(pIndex->GetIndexType());

#if defined(__GP_STRICT_VALIDATION)
	/// checkt that the index and 
	if(strcmp(GetZeroCurve()->GetCurrencyUnit()->GetCcyName(),
			pIndex->GetCurrencyUnit()->GetCcyName() ) != 0 )
	{
		CC_Ostringstream os;
		os	<< ARM_USERNAME << ": curve with " << GetZeroCurve()->GetCurrencyUnit()->GetCcyName()
			<< " ccy while index with " << pIndex->GetCurrencyUnit()->GetCcyName() << " ccy";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
	}
#endif

	/// create datestrip
	ARM_DateStrip* pFloatDateStrip = new ARM_DateStrip( startDate, endDate, floatResetFreq, floatDayCount, floatResetCal,
		fwdRule, intRule, stubRule, resetGap, floatPayFreq, floatPayGap, floatPayCal );

	delete floatResetCal;
	delete floatPayCal;

	return pFloatDateStrip;
}

//////////////////////////////////////////////////// ////////////////////////////////////////////				
///	Class  : ARM_MarkovSV
///	Routine: ComputeEquivalentCstShift
///	Returns: double Cst Shift
///	Action : This function is valid for 1 and 2 Factor Models
///			Returns the constant equivalent shift for the 
///			Shifted Heston Model
///
/// WARNING : TO BE UPDATED TO HAVE NBPOINTS PER YEAR!!!!
///			TO BE MERGED WITH CST VOL OF VOL CALCULATION
////////////////////////////////////////////////////////////////////////////////////////////////
double ARM_SVModels::ComputeEquivalentCstShift( 
		const string& curveName,
		double resetTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		double effetiveCstVolOfVol,
		double* globalVol) const
{
	double integratedVol = 0.;
	int n= GL_NBPOINTS_FIRST; // points for the first integration
	int n_bis = GL_NBPOINTS_SECOND; // points between the first points, for the calculation of ht

	// Constant coeffs for the volatility diffusion
	double volMrs		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetValueAtPoint(0);
	double volvol		= effetiveCstVolOfVol;

	double y_start,y_end,temp;
	double Sum = 0.0;
	double result = 0.0;
	double Sum_in = 0.0;
	double GlobalWeight = 0.0;    /// Integral (h�* sigma� , 0 , T)
	
	GaussLegendre_Coefficients c2(n);
	GaussLegendre_Coefficients c2_in(n_bis);

	///////////////////////Approach 2 ///////////////////
	y_start = 0.0;
	double weight,weight_in,y_end_in,t_in,Volt_in;

	for(int l=0; l<itsStoredData.itsNewSchedule.size();++l)
	{
		y_end = itsStoredData.itsNewSchedule[l]/K_YEAR_LEN;
		GaussLegendre_Coefficients coeff_out(&c2,y_start,y_end); 
		for(int k=0;k<n;++k)
		{
			Sum_in = 0.0;
			y_end_in = coeff_out.get_point(k);
			double time = y_end_in * K_YEAR_LEN;
			double Volt = ComputeEquivalentVolatilityt(curveName,time,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods);
			double Betat = ComputeEquivalentShiftt(curveName,time,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,Volt);
			double VolOfVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(time);


			GaussLegendre_Coefficients coeff_in(&c2_in,y_start,y_end_in); 


			double a = exp(-volMrs*y_end_in)*VolOfVol_t*VolOfVol_t;

			for(int m=0;m<n_bis;++m)
			{
				t_in = coeff_in.get_point(m);
				Volt_in = ComputeEquivalentVolatilityt(curveName,t_in * K_YEAR_LEN,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods);
				weight_in = coeff_in.get_weight(m);
				Sum_in += weight_in * (sinh(t_in*volMrs)*a/volMrs + 1.0) * Volt_in *Volt_in ;
			}
			Sum_in+= itsStoredData.itsBetaPart_1[l] + itsStoredData.itsBetaPart_2[l] * a;
			weight = coeff_out.get_weight(k);
			temp = weight * Volt * Volt;
			integratedVol+= temp;
			temp *= Sum_in;
			result+= Betat * temp;
			GlobalWeight+= temp;
		}
		y_start  = y_end;
	}
	*globalVol = integratedVol;

	return (result/GlobalWeight);
}


//////////////////////////////////////////////////// ////////////////////////////////////////////				
///	Class  : ARM_MarkovSV
///	Routine: ComputeEquivalentCstVolOfVol
///	Returns: double Cst Shift
///	Action : This function is valid for 1 and 2 Factor Models
///			Returns the constant equivalent shift for the 
///			Shifted Heston Model
///
/// WARNONG : TO BE UPDATED TO HAVE NBPOINTS PER YEAR!!!!
////////////////////////////////////////////////////////////////////////////////////////////////
double ARM_SVModels::ComputeEquivalentCstVolOfVol( 
		const string& curveName,
		double resetTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods) const
{
	int n= GL_NBPOINTS_FIRST; // points for the first integration
	int n_bis = GL_NBPOINTS_SECOND; // points between the first points, for the calculation of ht

	// Constant coeffs for the volatility diffusion
	double volMrs		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetValueAtPoint(0);

	double y_start,y_end;
	double Sum_in = 0.0;
	double result = 0.0;
	double GlobalWeight = 0.0;    /// Integral (h�* sigma� , 0 , T)
	
	GaussLegendre_Coefficients c2(n);
	GaussLegendre_Coefficients c2_in(n_bis);

	///////////////////////Approach 2 ///////////////////
	y_start = 0.0;
	double weight,weight_in,y_end_in,t_in,Volt_in;
	double Volt,temp,volvol_t,epsPart1,epsPart2;

	for(int l=0; l<itsStoredData.itsNewSchedule.size();++l)
	{
		y_end = itsStoredData.itsNewSchedule[l]/K_YEAR_LEN;		
		GaussLegendre_Coefficients coeff_out(&c2,y_start,y_end); 
		for(int k=0;k<n;++k)
		{
			Sum_in = 0.;
			y_end_in = coeff_out.get_point(k);
			Volt = ComputeEquivalentVolatilityt(curveName,y_end_in * K_YEAR_LEN,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods);
			temp = Volt * Volt * exp (- volMrs * y_end_in);
			GaussLegendre_Coefficients coeff_in(&c2_in,y_end_in,y_end); 
			for(int m=0;m<n_bis;++m)
			{
				t_in = coeff_in.get_point(m);
				Volt_in = ComputeEquivalentVolatilityt(curveName,t_in * K_YEAR_LEN,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods);
				weight_in = coeff_in.get_weight(m);
				Sum_in += weight_in * Volt_in * Volt_in * exp (-volMrs * t_in) ;
			}
			Sum_in+= itsStoredData.itsEpsPart_3[l+1];
			
			volvol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(y_end_in * K_YEAR_LEN);
			epsPart1 = computeEpsPart (y_start, y_end_in, volMrs, volvol_t);
			epsPart2 = computeEpsPart (y_start, y_end_in, volMrs, 1.0);
			
			epsPart1+= itsStoredData.itsEpsPart_1[l];
			epsPart2+= itsStoredData.itsEpsPart_2[l];

			weight = coeff_out.get_weight(k);
			double temp_in = weight  * temp * Sum_in;
			result+= temp_in * epsPart1;
			GlobalWeight+= temp_in * epsPart2;
		}
		y_start  = y_end;
	}
	return sqrt(result /GlobalWeight);
}

////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: ComputeInitialPoint
///	Returns: double 
///	Action : returns -g''(ksi)/g'(ksi) = 1/(2*ksi) + CstShift�/8
///			 ksi = V0 Integral ( sigma� , 0 , T)
////////////////////////////////////////////////////
double ARM_SVModels::ComputeInitialPoint( 
	double cstShift)
{
	
	double ksi = itsIntegratedVol;

#if defined(__GP_STRICT_VALIDATION)
	if(ksi<K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Integrated Vol of Vol is Null!");
#endif
	return (1/(2* ksi) + cstShift*cstShift/8);
}
 

///////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: ComputePhi
///	Returns: double 
///	Action : returns the Value of Phi at the InitialPoint Previously calculated
///			We a RungeKutta to solve the 2D Ordinary Differential Equation
///////////////////////////////////////////////////////////////////////////////////
double ARM_SVModels::ComputePhi(const ARM_ODEFunc& test,
							double cstShiftValue,
							double resetTime) const

{
	double V0			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol))).GetValueAtPoint(0);// Equivalent Vol Of Vol	
	double A0T = 0.;
	double B0T = 0.;

	// We use Runge Kutta
	int nvar = 2; // 2 Variables for AtT and BtT
	double eps = ODE_PRECISION;
	double tiny = 1.0e-30;

	vector<std::vector<double>> solverVars(RK5_NBTMP_VAR > RK4_NBTMP_VAR ? RK5_NBTMP_VAR : RK4_NBTMP_VAR);
	for(size_t i=0;i<solverVars.size();++i)
		solverVars[i].resize(nvar);

	std::vector<double> yInOut(nvar,0.0); // start at 0

	double x1 = resetTime/K_YEAR_LEN; // Start
	double x2 = 0;// End Time

	double h1 = - 0.3;
	double hmin = 0.0;
	int nok = 0;
	int nbad = 0;
	int kmax = 1000;
	int kount=0;

	std::vector<double>& xp = new std::vector<double>(kmax * nvar,0.0);

	double dxsav =0.3;
	int nbSteps;

	odeint(&yInOut, nvar, x1, x2, eps,tiny,h1,hmin,kmax,dxsav,test,solverVars,&nok,&nbad,&kount,nbSteps,xp,NULL);

	A0T = yInOut[0];
	B0T = yInOut[1];

	return exp(A0T -  B0T);
}

///////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: ComputePhi_Zero
///	Returns: double 
///	Action : returns the Value of Phi_Zero at the InitialPoint Previously calculated
///			We a The Standard Riccati Solver for Piece-Wise Coefficients
//////////////////////////////////////////////////////////////////////////////////
double ARM_SVModels::ComputePhi_Zero( 
	double cstShift,
	double resetTime) const
{
	double V0 = ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol))).GetValueAtPoint(0);// Equivalent Vol Of Vol	
	double A0T = 0.;
	double B0T = 0.;


	double T = resetTime/K_YEAR_LEN;

	ARM_RiccatiMSV test((*GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam().GetCurve()),
						(*GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()),
						(*GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()),
						(*GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()));


	A0T = test.A(0,T);
	B0T = test.B0T(T);
	return exp(A0T - V0 * B0T);
}

///////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: ComputeMu
///	Returns: double 
///	Action : returns the Value of Mu (input for the Riccati)
///////////////////////////////////////////////////////////////////////////////////
double ARM_SVModels::ComputeMu(double cstShiftValue,
							double integratedVol) const
{
	return (1/(2*integratedVol)+cstShiftValue*cstShiftValue/8);
}



///////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: UpdateRiccatiAt
///	Returns: double 
///	Action : returns the Value of Ct at time t for the Riccati (non Analytical) Equation
///////////////////////////////////////////////////////////////////////////////////
double ARM_SVModels::UpdateRiccatiCt(const string& curveName,
									double t,
									double floatStartTime,
									double floatEndTime,
									const std::vector<double>& fixPayTimes,
									const std::vector<double>& fixPayPeriods)
{
	double result = ComputeEquivalentVolatilityt(curveName,t * K_YEAR_LEN,floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods);
	return result;
}
double ARM_SVModels::UpdateRiccatiBt(double t)
{
	double result = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(t * K_YEAR_LEN);
	return result;
}
double ARM_SVModels::UpdateRiccatiAt(double t)
{
	double result = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(t * K_YEAR_LEN);
	result = 0.5 *result * result;
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: ComputeEquivalentCstVolatility
///	Returns: double Cst Vol
///	Action : This function is valid for 1 and 2 Factor Models
///			Returns the constant equivalent volatility for the 
///			Shifted Heston Model
////////////////////////////////////////////////////
double ARM_SVModels::ComputeEquivalentCstVolatility( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_PricingStatesPtr& states) const
{
	return 0.0;
}

////////////////////////////////////////////////////
///	Class   : ARM_SVModels
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : void
///	Action  : VolatilitiesAndCorrelationTimesSteps for PDE
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_SVModels::VolatilitiesAndCorrelationTimesSteps() const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VolatilitiesAndCorrelationTimesSteps : unimplemented function for ARM_MarkovSV Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_SVModels
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : Default implementation is no markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_SVModels::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix(numMethodStates->rows(), numMethodStates->cols(), 0.0 ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SVModels::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in MSV1F model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the MSV1F model!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_SVModels
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SVModels::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
   	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Libor : unimplemented function for ARM_MarkovSV Model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SVModels::VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaCaplet : unimplemented function for ARM_MarkovSV Model!");
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SVModels::VanillaDigital(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
        double fwdPeriod,
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaDigital : unimplemented function for ARM_MarkovSV Model!");
}



////////////////////////////////////////////////////
///	Class  : ARM_SVModels
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SVModels::VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const
{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSwaption : unimplemented function for ARM_MarkovSV Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV
///	Routine : ComputeModelTimes
///	Returns : an empty vector since in HW there is not
///				such a thing as model times
////////////////////////////////////////////////////
std::vector<double>& ARM_SVModels::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// since there is no concept of model time
	/// returns an empty vector
	return new std::vector<double>(0);
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SVModels::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"MCModelStatesFromToNextTime : unimplemented function for ARM_MarkovSV Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SVModels::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_MarkovSV Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_SVModels::VanillaSpreadOptionLet(const string& curveName,
													double evalTime,
													int callPut,
													double startTime,
													double endTime,
													double resetTime,
													double payTime,
													double payPeriod,
													double notional,
													double coeffLong,
													double coeffShort,
													const std::vector<double>& strikes,
													double swapLongFloatStartTime,
													double swapLongFloatEndTime,
													const std::vector<double>& swapLongFixPayTimes,
													const std::vector<double>& swapLongFixPayPeriods,
													double swapShortFloatStartTime,
													double swapShortFloatEndTime,
													const std::vector<double>& swapShortFixPayTimes,
													const std::vector<double>& swapShortFixPayPeriods,
													const ARM_PricingStatesPtr& states) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSpreadOption : unimplemented function for ARM_MarkovSV Model!");
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

