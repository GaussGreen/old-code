#include "stdafx.h"
#include "Heston.h"






heston_fitter::heston_fitter( MlEqAssetHandle hudly, long dateMaturity, int nPoints ):
CStaticReplication(hudly, dateMaturity, nPoints),
m_hUnderlying(hudly)
{
	m_nToday = m_hUnderlying->GetDateHandle()->GetDate();

	m_dateShort	 = MlEqMaths::Min(m_nToday + 90, (m_dateMaturity+m_nToday)/2) ;	// 1m
	m_dateMiddle = m_dateMaturity;
	m_dateLong	 = MlEqMaths::Max(m_nToday + 2000, 2*m_dateMaturity - m_nToday); // 10y

	m_params.v0		 = 0.;
	m_params.vas	 = 0.;
	m_params.correl	 = 0.;
	m_params.kappa	 = 0.;
	m_params.vovol	 = 0.;
}


heston_fitter::heston_fitter( MlEqAssetHandle hudly,
							 long dateShort,
							 long dateMaturity,
							 long dateLong,
							 int nPoints ):
CStaticReplication(hudly, dateMaturity, nPoints),
m_hUnderlying(hudly)
{
	m_nToday = m_hUnderlying->GetDateHandle()->GetDate();

	m_dateShort	 = dateShort;
	m_dateMiddle = dateMaturity;
	m_dateLong	 = dateLong;

	m_params.v0		 = 0.;
	m_params.vas	 = 0.;
	m_params.correl	 = 0.;
	m_params.kappa	 = 0.;
	m_params.vovol	 = 0.;
}


void heston_fitter::calibrate( GVector<MlEqStrikeHandle> vStrikes )
{
	m_params.vovol = 0.5;		// first guess for calibration
	m_params.correl = -0.7;

	bool stop = false;

	while( !stop )
	{
		calibrateToVarianceSwap();
		stop =  calibrateToVanilla( vStrikes ) ;
		checkMaturities( stop );
	}
}


void heston_fitter::checkMaturities(bool& stop)
{
	if( stop )	return;

	if( m_dateShort == (m_dateMaturity+m_nToday)/2 )
	{
		stop = true;	// last case
	}
	else
	{
		m_dateShort = (m_dateMaturity+m_nToday)/2 ;	
		stop = false;
	}
}


void  heston_fitter::ObjectiveFcn(double* gVals,double* xVals)
{
	if( m_calibVanilla )
		ObjectiveFcnVanilla(gVals, xVals);
	else
		ObjectiveFcnVS(gVals, xVals);
}


void heston_fitter::initCalibrationInputs(GVector<MlEqStrikeHandle> vStrikes)
{
	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	
	m_mat	= hDate->GetYearFraction(m_dateMaturity);
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	m_forward =   m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);

	MlEqVolatilityStructureHandle	hVolatilityStructure = m_hUnderlying->GetVolatilityStructure();
	double spot = m_hUnderlying->GetSpot(nToday);	
	
	m_dataSize = vStrikes.getsize();
	if( m_dataSize )
	{
		m_normStrikes.resize(0);
		m_strikes.resize(m_dataSize,0.);
		MlEqStrike stk;	

		for (int i = 0; i < m_dataSize; i++)
		{
			MlEqStrike::convertStrikes(stk,*(vStrikes[i]) );
			m_strikes[i] = stk.m_strike;
		}
	}
	else
	{
		m_dataSize = 5;
		m_normStrikes.resize(m_dataSize,0.);
		m_strikes.resize(m_dataSize,0.);
		m_normStrikes[0] = -1.0;
		m_normStrikes[1] = -0.5;
		m_normStrikes[2] =  0.0;
		m_normStrikes[3] =  0.5;
		m_normStrikes[4] =  1.0;

		double vol0 = hVolatilityStructure->getFutureVol(MlEqStrike(m_forward), nToday, m_dateMaturity, spot);

		for(int i=0; i<5; i++){
			m_strikes[i] = m_forward* exp(vol0*sqrt(m_mat)* m_normStrikes[i]) ;	
		}

	}
	
	m_prices.resize(m_dataSize,0.);
	m_bsVega.resize(m_dataSize,0.);

	for(int i=0; i<m_dataSize; i++)
	{
		double strike = m_strikes[i];
		double vol  = hVolatilityStructure->getFutureVol(MlEqStrike(strike), nToday, m_dateMaturity, spot);
	
		m_prices[i]	=	 ::Bs(m_forward, vol, m_mat, strike, 1.0, 1.);	
		m_bsVega[i] =	(::Bs(m_forward, vol+0.01, m_mat, strike, 1.0, 1.) - m_prices[i] ) * 100. ;
	}

}



bool heston_fitter::calibrateToVanilla( GVector<MlEqStrikeHandle> vStrikes )
{
	m_calibVanilla	=	true;
	m_pricer = new heston_pricer(m_params);
	initCalibrationInputs(vStrikes);



	CMatrix NonZeroPartialDerivativesSpecification;

	double FinalTolerance		= 1e-8;
	double StoppingTolerance	= 1e-8;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	int returnCode;

	iVector linearVariableIndex(2,0);
	CMatrix ObjectiveBounds;	
	CVector initialGuess(2);
	CMatrix xValBounds(2,2);
	CVector x;

	initialGuess[0]	= m_pricer->m_cov ;	
	initialGuess[1]	= m_pricer->m_var ;	
	
	xValBounds[0][0] = -3.;
	xValBounds[0][1] =  3.;
	xValBounds[1][0] =  0.01;
	xValBounds[1][1] =  3.;

	initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	x = initialGuess ;	
	solve(x, 0, returnCode);

	m_pricer->m_cov	 = x[0];
	m_pricer->m_var	 = x[1];

	double err = checkFcn();

///	double check;
//	ObjectiveFcn(&check, x);

	m_params.vovol = sqrt(m_pricer->m_var);
	m_params.correl = m_pricer->m_cov/m_params.vovol ;

	return ( err < 0.0010 );	// 10 bp on implied vol	
}



void heston_fitter::ObjectiveFcnVanilla(double* gVals,double* xVals)
{
	*gVals = 0.0;
	
	m_pricer->m_cov =  xVals[0];
	m_pricer->m_var =  xVals[1];

	for(int i=0; i<m_dataSize; i++)
	{
		double price = m_pricer->call_price(m_mat, m_forward, m_strikes[i]);
		price -= m_prices[i] ;
		price /= m_bsVega[i] ;
		*gVals += fabs(price);			
	}
}


double heston_fitter::checkFcn()
{
	double err = 0.0;

	for(int i=0; i<m_dataSize; i++)
	{
		double price = m_pricer->call_price(m_mat, m_forward, m_strikes[i]);
		price -= m_prices[i] ;
		price /= m_bsVega[i] ; 

		err = std::max( err, fabs(price) );
	}
	return err;
}




void heston_fitter::calibrateToVarianceSwap()
{
	m_calibVanilla	=	false;

	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	

	m_t.resize(3,0.);
	m_t[0]  = hDate->GetYearFraction(m_dateShort);
	m_t[1]  = hDate->GetYearFraction(m_dateMiddle);	
	m_t[2]  = hDate->GetYearFraction(m_dateLong);

	m_vs.resize(3,0.);

	m_dateMaturity	=	m_dateShort;
	m_mat			=	m_t[0];
	m_vs[0]			=	CStaticReplication::fComputePrice() / m_t[0] ;

	m_dateMaturity   =	m_dateMiddle;
	m_mat			=	m_t[1];
	m_vs[1]			=	CStaticReplication::fComputePrice() / m_t[1] ;

	m_dateMaturity	=	m_dateLong;
	m_mat			=	m_t[2];
	m_vs[2]			=	CStaticReplication::fComputePrice() / m_t[2] ;


	m_dateMaturity	=	m_dateMiddle ;


	CMatrix NonZeroPartialDerivativesSpecification;
	double FinalTolerance		= 1e-12;
	double StoppingTolerance	= 1e-12;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	
	CVector initialGuess(3);
	initialGuess[0]	= m_vs[2] * m_t[2] ;
	initialGuess[1]	= (m_vs[0] * m_t[0] - m_vs[2] * m_t[2])/0.4;
	initialGuess[2]	= 0.4;
	
	CMatrix xValBounds(3,2);
	xValBounds[0][0] = 0.01;
	xValBounds[0][1] = 1.  ;
	xValBounds[1][0] = -5.0;
	xValBounds[1][1] = 5.0 ;
	xValBounds[2][0] = 0.0 ;
	xValBounds[2][1] = 10. ;

	iVector linearVariableIndex(3,0);
	CMatrix ObjectiveBounds;
	
	LSGRGSolver::initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	CVector x = initialGuess ;
	int returnCode;;
	LSGRGSolver::solve(x, 0, returnCode);

	m_params.kappa	 = x[2];
	m_params.vas	 = x[0];
	m_params.v0		 = x[2]*x[1] + x[0];
}



void heston_fitter::ObjectiveFcnVS(double* gVals,double* xVals)
{
	*gVals = 0.0;	

	for(int i=0; i<3; i++)
	{
		double vs_h = xVals[0] + xVals[1] * (1.-exp(-xVals[2]*m_t[i]) )/m_t[i] ;
		vs_h -= m_vs[i] ;
		vs_h *= vs_h ;
		*gVals += vs_h;
	}
}

//////////////////////////////////////////////////



heston_jump_fitter::heston_jump_fitter( MlEqAssetHandle hudly, long dateMaturity, int nPoints ):
heston_fitter(hudly, dateMaturity, nPoints)
{
	m_lambda	=	m_jump	= 0.;
}


heston_jump_fitter::heston_jump_fitter( MlEqAssetHandle hudly,
										 long dateShort,
										 long dateMaturity,
										 long dateLong,
										 int nPoints ):
heston_fitter(hudly, dateShort, dateMaturity, dateLong, nPoints)
{
	m_lambda	=	m_jump	= 0.;
}


void heston_jump_fitter::calibrate( GVector<MlEqStrikeHandle> vStrikes )
{
	m_params.vovol = 0.001;		// first guess for calibration
	m_params.correl = 0.;

	calibrateJumpsToSkew(vStrikes);
//	calibrateToVarianceSwap();
	calibrateToVanilla(vStrikes);
}

void heston_jump_fitter::initCalibrationInputs(GVector<MlEqStrikeHandle> vStrikes)
{
//	heston_fitter::initCalibrationInputs(vStrikes);
//	return;

	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	
	m_mat	= hDate->GetYearFraction(m_dateMaturity);
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	m_forward =   m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);

	MlEqVolatilityStructureHandle	hVolatilityStructure = m_hUnderlying->GetVolatilityStructure();
	double spot = m_hUnderlying->GetSpot(nToday);	
	
	m_dataSize = vStrikes.getsize();
	if( m_dataSize )
	{
		m_normStrikes.resize(0);
		m_strikes.resize(m_dataSize,0.);
		MlEqStrike stk;	

		for (int i = 0; i < m_dataSize; i++)
		{
			MlEqStrike::convertStrikes(stk,*(vStrikes[i]) );
			m_strikes[i] = stk.m_strike;
		}
	}
	else
	{
		double vol0 = hVolatilityStructure->getFutureVol(MlEqStrike(m_forward), nToday, m_dateMaturity, spot);

		if(m_calibJumps)
		{
			m_dataSize = 4;
			m_normStrikes.resize(m_dataSize,0.);
			m_strikes.resize(m_dataSize,0.);
			m_normStrikes[0] = -2.5;
			m_normStrikes[1] = -1.5;
			m_normStrikes[2] = -0.5;
			m_normStrikes[3] =  0.5;

			m_bsVol = vol0;

		}
		else
		{
			m_dataSize = 7;
			m_normStrikes.resize(m_dataSize,0.);
			m_strikes.resize(m_dataSize,0.);
			m_normStrikes[0] = -1.5;
			m_normStrikes[1] = -0.5;
			m_normStrikes[2] =  0.0;
			m_normStrikes[3] =  0.5;
			m_normStrikes[4] =  1.0;

			m_normStrikes[5] =  -2.0;
			m_normStrikes[6] =  -2.5;
		}

	
		for(int i=0; i<m_dataSize; i++){
			m_strikes[i] = m_forward* exp(vol0*sqrt(m_mat)* m_normStrikes[i]) ;	
		}

	}
	
	m_vols.resize(m_dataSize,0.);
	m_prices.resize(m_dataSize,0.);
	m_bsVega.resize(m_dataSize,0.);

	for(int i=0; i<m_dataSize; i++)
	{
		double strike = m_strikes[i];
		double vol  = hVolatilityStructure->getFutureVol(MlEqStrike(strike), nToday, m_dateMaturity, spot);
	
		m_vols[i]	=	vol ;
		m_prices[i]	=	 ::Bs(m_forward, vol, m_mat, strike, 1.0, 1.);	
		m_bsVega[i] =	(::Bs(m_forward, vol+0.01, m_mat, strike, 1.0, 1.) - m_prices[i] ) * 100. ;
	}

}



bool heston_jump_fitter::calibrateToVanilla( GVector<MlEqStrikeHandle> vStrikes )
{
	m_calibJumps	=	false;
	m_calibVanilla	=	true;
	m_jump_pricer = new heston_jump_pricer(m_params, m_lambda, m_jump);
	initCalibrationInputs(vStrikes);



	CMatrix NonZeroPartialDerivativesSpecification;

	double FinalTolerance		= 1e-8;
	double StoppingTolerance	= 1e-8;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	int returnCode;

	iVector linearVariableIndex(4,0);
	CMatrix ObjectiveBounds;	
	CVector initialGuess(4);
	CMatrix xValBounds(4,2);
	CVector x;

	initialGuess[0]	= -0.25*0.25 ;	
	initialGuess[1]	= 0.25*0.25;
	initialGuess[2]	= 0.1;
	initialGuess[3]	= m_bsVol*m_bsVol ;

	xValBounds[0][0] = -3.;
	xValBounds[0][1] =  3.;
	xValBounds[1][0] =  0.01;
	xValBounds[1][1] =  3.;
	xValBounds[2][0] =  0.;
	xValBounds[2][1] =  10.;
	xValBounds[3][0] =  0.0001;
	xValBounds[3][1] =  1.;

	initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	x = initialGuess ;	
	solve(x, 0, returnCode);

	m_jump_pricer->m_cov	 = x[0];
	m_jump_pricer->m_var	 = x[1];
	m_jump_pricer->m_lambda	 = x[2];
	m_jump_pricer->m_vas	 = m_jump_pricer->m_v0		 = x[3];
	m_jump_pricer->m_kappa = 0.;

	double err ;
	ObjectiveFcn(&err, x);
	err /= double(m_dataSize) ;

	m_params.vovol = sqrt(m_jump_pricer->m_var);
	m_params.correl = m_jump_pricer->m_cov/m_params.vovol ;
	m_lambda = 	m_jump_pricer->m_lambda ;
	m_jump	 =  m_jump_pricer->m_jump ;
	m_params.v0		=		m_jump_pricer->m_v0;
	m_params.vas	=		m_jump_pricer->m_vas;
	m_params.kappa	=		m_jump_pricer->m_kappa;

	return ( err < 0.0010 );	// 10 bp on implied vol	
}




void heston_jump_fitter::calibrateJumpsToSkew(GVector<MlEqStrikeHandle> vStrikes)
{
	m_calibJumps	=	true;
	m_calibVanilla	=	true;

	initCalibrationInputs(vStrikes);



	CMatrix NonZeroPartialDerivativesSpecification;

	double FinalTolerance		= 1e-8;
	double StoppingTolerance	= 1e-8;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	int returnCode;

	iVector linearVariableIndex(2,0);
	CMatrix ObjectiveBounds;	
	CVector initialGuess(2);
	CMatrix xValBounds(2,2);
	CVector x;

	initialGuess[0]	= 0.1;//m_lambda ;	
	initialGuess[1]	= -0.2;//m_jump ;	

	xValBounds[0][0] =  0.001;
	xValBounds[0][1] =  10.;
	xValBounds[1][0] =  -1.;
	xValBounds[1][1] =  2.;

	initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	x = initialGuess ;	
	solve(x, 0, returnCode);
}


void heston_jump_fitter::ObjectiveFcnVanilla(double* gVals,double* xVals)
{
	*gVals = 0.0;

	if(m_calibJumps)
	{
		m_lambda =  xVals[0];
		m_jump	 =  xVals[1];

		double atmvol = m_bsVol;

		double price = MertonPrice(m_mat, m_forward, m_forward);
		atmvol -= (price - m_prices[2])/m_bsVega[2];

		for(int i=0; i<m_dataSize; i++)
		{
			double price = MertonPrice(m_mat, m_forward, m_strikes[i]);
			price -= m_prices[i] ;
			price /= m_bsVega[i] ;
			
			double dslope = price/atmvol - (atmvol/m_bsVol - 1.)*m_vols[i]/m_bsVol ;

			*gVals += fabs(dslope);			
		}
	}
	else
	{
	
		m_jump_pricer->m_cov	 =  xVals[0];
		m_jump_pricer->m_var	 =  xVals[1];
		m_jump_pricer->m_lambda	 =  xVals[2];

		m_jump_pricer->m_v0		 =  xVals[3];
		m_jump_pricer->m_vas	 =  xVals[3];
		m_jump_pricer->m_kappa	 =  0.0;//xVals[6];

		for(int i=0; i<m_dataSize; i++)
		{
			double price = m_jump_pricer->call_price(m_mat, m_forward, m_strikes[i]);
			price -= m_prices[i] ;
			price /= m_bsVega[i] ;
			*gVals += fabs(price);			
		}
	}
}


double heston_jump_fitter::checkFcn()
{
	double err = 0.0;

	for(int i=0; i<m_dataSize; i++)
	{
		double price = m_jump_pricer->call_price(m_mat, m_forward, m_strikes[i]);
		price -= m_prices[i] ;
		price /= m_bsVega[i] ; 

		err = std::max( err, fabs(price) );
	}
	return err;
}

double heston_jump_fitter::MertonPrice(double mat, double forward, double strike)
{
	double lambdat = m_lambda*m_mat;
	double drift	= -lambdat*(exp(m_jump)-1.);
	int factorial = 1;
	double price = 0.;

	for(int njump = 0; njump<15; njump++)
	{
		double bsPrice = ::Bs(forward*exp(njump*m_jump+drift), m_bsVol, mat, strike, 1., 1., 0.);
		double dprice = bsPrice * pow(lambdat, njump) / factorial;
		price += dprice;

		if(dprice/forward < 1e-6 )	break;

		if(njump>0)
			factorial *=  njump;
	}

	return price*exp(-lambdat) ;
}




void heston_jump_fitter::calibrateToVarianceSwap()
{
	m_calibVanilla	=	false;

	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	

	m_t.resize(4,0.);
	m_vs.resize(4,0.);

//	m_t[1]  = hDate->GetYearFraction(m_dateMiddle);	
//	m_t[2]  = hDate->GetYearFraction(m_dateLong);

	
	m_dateMaturity	=	(m_dateMiddle + 2*m_nToday)/3;
	m_t[0]			=	hDate->GetYearFraction(m_dateMaturity);
	m_mat			=	m_t[0];
	m_vs[0]			=	CStaticReplication::fComputePrice() / m_t[0] ;

	m_dateMaturity	=	(2*m_dateMiddle + m_nToday)/3;
	m_t[1]			=	hDate->GetYearFraction(m_dateMaturity);
	m_mat			=	m_t[1];
	m_vs[1]			=	CStaticReplication::fComputePrice() / m_t[1] ;

	m_dateMaturity	=	m_dateMiddle;
	m_t[2]			=	hDate->GetYearFraction(m_dateMaturity);
	m_mat			=	m_t[2];
	m_vs[2]			=	CStaticReplication::fComputePrice() / m_t[2] ;

	m_dateMaturity	=	m_dateLong;
	m_t[3]			=	hDate->GetYearFraction(m_dateMaturity);
	m_mat			=	m_t[3];
	m_vs[3]			=	CStaticReplication::fComputePrice() / m_t[3] ;


//	double vsRef = m_vs[1];
//	for(int i=0; i<3; i++)
//		m_vs[i] /= vsRef ;


	m_dateMaturity	=	m_dateMiddle ;


	CMatrix NonZeroPartialDerivativesSpecification;
	double FinalTolerance		= 1e-12;
	double StoppingTolerance	= 1e-12;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	
	CVector initialGuess(4);
	initialGuess[0]	= m_vs[3] * m_t[3] ;
	initialGuess[1]	= (m_vs[0] * m_t[0] - m_vs[3] * m_t[3])/0.4;
	initialGuess[2]	= 0.4;
	initialGuess[3]	= 0.1 ;
	
	CMatrix xValBounds(4,2);
	xValBounds[0][0] = 0.01;
	xValBounds[0][1] = 1.  ;
	xValBounds[1][0] = -5.0;
	xValBounds[1][1] = 5.0 ;
	xValBounds[2][0] = 0.0 ;
	xValBounds[2][1] = 10. ;
	xValBounds[3][0] = 0.0 ;
	xValBounds[3][1] = 10. ;


	iVector linearVariableIndex(4,0);
	CMatrix ObjectiveBounds;
	
	LSGRGSolver::initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	CVector x = initialGuess ;
	int returnCode;;
	LSGRGSolver::solve(x, 0, returnCode);

	m_params.kappa	 = x[2];
	m_params.vas	 = x[0];
	m_params.v0		 = x[2]*x[1] + x[0];
	m_lambda		 = x[3];
}


void heston_jump_fitter::ObjectiveFcnVS(double* gVals,double* xVals)
{
	*gVals = 0.0;	

	for(int i=0; i<4; i++)
	{
		double vs_h = xVals[0] + xVals[1] * (1.-exp(-xVals[2]*m_t[i]) )/m_t[i] ;
		vs_h +=	2.*xVals[3]*(exp(m_jump)-1.-m_jump);
		vs_h -= m_vs[i] ;
		vs_h *= vs_h ;
		*gVals += vs_h;
	}
}


//////////////////////////////////////////////////


heston_pricer::heston_pricer(HestonParameters params)
{
	m_mat	= 0.;

	m_kappa = params.kappa ;
	m_v0	= params.v0 ;
	m_vas	= params.vas;

	m_var	= params.vovol;
	m_cov	= m_var * params.correl ;
	m_var  *= m_var;

	m_nPoints		= 30;
	m_gaussWeights	= CVector(m_nPoints, 0.0);
	m_gaussPoints	= CVector(m_nPoints, 0.0);

	MlEqMaths::dGauleg(0.0, 1., m_gaussPoints, m_gaussWeights, m_nPoints, false);
}


double heston_pricer::call_price(double mat, double forward, double strike)
{
	m_mat = mat;
	double logstrike = log(strike/forward);

	double price = 0.0;
	double P1 = 0.0, P2 = 0.0 ;
	complex i(0., 1.);

	double X = searchUpperBound();


	for ( int j = 0 ; j < m_nPoints; j++ )
	{
		double omega = X * m_gaussPoints[j];		
		complex o = exp(-(omega*logstrike)*i);

		complex temp = log_fft(omega);
		double dP1 = (-i*o*temp/omega).x ;

		temp = log_fft(omega-i)  ;
		double dP2 = (-i*o*temp/omega).x ;


		P1 += m_gaussWeights[j]* dP1 * X;
		P2 += m_gaussWeights[j]* dP2 * X;
	}
	

	price =  forward * ( 0.5 + P2 / 3.141592653590 );
	price -= strike  * ( 0.5 + P1 / 3.141592653590 );

	return price ;
}

double heston_pricer::searchUpperBound()
{
	complex temp;
	complex i(0., 1.);
	double ud = 1.;
	double bound = 100.;
	double dX = 50.;

	double dXmax = 10.;

	while( dX > dXmax )
	{

		double omega = bound;	

		temp = log_fft(omega)/omega;
		double eps = norm(temp) ;
		temp = log_fft(omega-i) ;
		eps = std::max( eps, norm(temp) ) ;

		ud = (eps < 1e-15)? -1.:1. ;
		bound += ud*dX;

		dX *= 0.5;
		
		if( bound < 10. )		dXmax = 0.4;
		else if( bound < 40. )	dXmax = 1.;
		else if( bound < 80. )	dXmax = 6.;
	}

	bound += 2.*dX;

	return bound;
}


complex heston_pricer::log_fft(complex c)
{
	complex i(0.,1.);

	complex temp = -m_kappa + i*c*m_cov;
	complex d = temp*temp;

	temp = c*(c + i);

	d = sqrt( d + temp*m_var );

	complex g =  - d - i*m_cov*c  + m_kappa ;
	temp = g + 2.*d ;
	g = temp / g;

	complex C,D;

	if( d.x*m_mat > 50.)
	{
		D = temp / (g*m_var);
		C = (temp - 2.*d)*m_mat - 2.*log(g/(g-1.));
	}
	else
	{
		complex temp2 =  exp(d*m_mat);
		complex temp3 =  1. - g*temp2 ;
		temp2 = 1. - temp2;

		D = (temp*temp2) / (temp3*m_var);

		temp3 = temp3 / (1.-g);
		temp3 = log( temp3 );

		C = (temp*m_mat - 2.*temp3) * m_kappa*m_vas/m_var ;
	}

	temp = C + D*m_v0;


	if( temp.x < -60.)	return 0.0;

	temp = exp( temp ) ;

	return temp;
}



double heston_pricer::implied_vol(double fwd, double mat, double strike)
{
	double price = call_price(mat, fwd, strike);
	
	double accuracy   = 1e-16;
	double lower_x    = 0.01; 
	double upper_x    = 1.5; 

	int rootfind_flag = eROOTFIND_BRENT_NOGROW;

	double vol =	MlEqBSImpliedVol(
								price,
								fwd,
								mat,
								strike,
								1,1,rootfind_flag ,	 
								accuracy   ,		 
								lower_x    ,		 
								upper_x);	

	return vol;
}


heston_jump_pricer::heston_jump_pricer(HestonParameters params, double lambda, double jump)
:heston_pricer(params),
m_lambda(lambda),
m_jump(jump)
{
}


complex heston_jump_pricer::log_fft(complex c)
{
	complex i(0.,1.);

	complex jump_fft = exp(i*c*m_jump) - 1. - i*c*(exp(m_jump) - 1.) ;
	jump_fft = exp(m_mat*m_lambda*jump_fft);

	complex heston_fft = heston_pricer::log_fft(c);

	return heston_fft * jump_fft;
}



/*
***
****
*****
****
***
*/
/*
***
****
*****
******
*******
********
*******
******
*****
****
***
**
*/




heston_term_pricer::heston_term_pricer(HestonParametersTermStructure params):
heston_pricer(params.hParams[0]),
m_t(params.t)
{
	m_nLastSlice = m_nSlice = m_t.getsize();

	m_vKappa.resize(m_nSlice,0.);
	m_vCov.resize(m_nSlice,0.);	
	m_vVar.resize(m_nSlice,0.);	

	m_dt.resize(m_nSlice, m_t[0]);

	for(int k=0; k<m_nSlice; k++)
	{
		m_vKappa[k] = params.hParams[k].kappa ;
		m_vCov[k]	= params.hParams[k].correl * params.hParams[k].vovol ;
		m_vVar[k]	= params.hParams[k].vovol  * params.hParams[k].vovol ;

		if(k>0)
			m_dt[k] = m_t[k] - m_t[k-1];
	}

	m_v0  = params.hParams[0].v0 ;
	m_vas = params.hParams[0].vas ;
}

double heston_term_pricer::call_price(double mat, double forward, double strike)
{
	updateTimeSlices(mat);
	double price = heston_pricer::call_price(mat, forward, strike);
	return price;
}

void heston_term_pricer::updateTimeSlices(double mat)
{
	m_mat = mat;

	for(int t=0; t<m_nSlice; t++)
	{
		if( mat < m_t[t]+1e-6 )
		{
			m_nLastSlice = t+1;

			if(t > 0 )				
				m_dt[t] = mat - m_t[t-1];
			else
				m_dt[0] = mat;

			return;
		}	
	}

	// case where mat is larger than all maturities

	m_nLastSlice = m_nSlice;
	if(m_nSlice > 1)
		m_dt[m_nSlice-1] = mat - m_t[m_nSlice-2];
	else
		m_dt[0] = mat;
}


complex heston_term_pricer::log_fft(complex c)
{
	complex i(0.,1.);	

	complex temp = -m_vKappa[0] + i*c*m_vCov[0];
	complex d = temp*temp;

	temp = c*(c + i);

	d = sqrt( d + temp*m_vVar[0] );

	complex g =  - d - i*m_vCov[0]*c  + m_vKappa[0] ;
	temp = g + 2.*d ;
	g = temp / g;

	complex Cm,Dm;
	d = int_d(c);

	if( d.x > 50.)
	{
		Dm = temp / (g*m_vVar[0]);
	}
	else
	{
		complex temp2 =  exp(d);
		complex temp3 =  1. - g*temp2 ;
		temp2 = 1. - temp2;

		Dm = (temp*temp2) / (temp3*m_vVar[0]);
	}

	Cm = C(c);

	temp = Cm + Dm*m_v0;


	if( temp.x < -60.)	return 0.0;

	temp = exp( temp ) ;

	return temp;
}

complex heston_term_pricer::int_d(complex c)
{
	complex sum(0.,0.);
	complex i(0.,1.);

	for(int n=0; n<m_nLastSlice; n++)
	{
		complex tmp = -m_vKappa[n] + i*c*m_vCov[n];
		complex d = tmp*tmp;
		tmp = c*(c + i);
		d = sqrt( d + tmp*m_vVar[n] );

		sum = sum + d*m_dt[n];
	}
	return sum;
}


complex heston_term_pricer::int_d(complex c, int start, int end)
{
	complex sum(0.,0.);
	complex i(0.,1.);

	end = std::min(end, m_nLastSlice);

	for(int n=start; n<end; n++)
	{
		complex tmp = -m_vKappa[n] + i*c*m_vCov[n];
		complex d = tmp*tmp;
		tmp = c*(c + i);
		d = sqrt( d + tmp*m_vVar[n] );

		sum = sum + d*m_dt[n];
	}
	return sum;
}




complex heston_term_pricer::C(complex c)
{
	complex sum(0.,0.);
	complex i(0.,1.);

	for(int n=0; n<m_nLastSlice; n++)
	{
		complex tmp = -m_vKappa[n] + i*c*m_vCov[n];
		complex d = tmp*tmp;

		d = sqrt( d + c*(c + i)*m_vVar[n] );

		complex g =  - d - tmp; 
		tmp = g + 2.*d ;
		g = tmp / g;
		double dt = m_dt[n];

		complex intd_term = int_d(c, n+1, m_nLastSlice);
		complex log_term(0.,0.);
	
		if( intd_term.x > 50.)
		{			
			log_term = dt * d;		
		}
		else if( dt * d.x > 50.)
		{
			log_term = d*dt - log( exp(intd_term)-1/g) ;
		}
		else
		{
			intd_term = g*exp( intd_term );
			log_term = log( (1-intd_term*exp(d*dt))/(1-intd_term) );
		}

		tmp = tmp*dt - 2.*log_term ;

		sum = sum + m_vKappa[n]*m_vas/m_vVar[n]*tmp;
	}
	return sum;
}







heston_term_fitter::heston_term_fitter( MlEqAssetHandle udly, GVector<long> dateMaturities ):
heston_fitter(udly, dateMaturities[0]),
m_dateMaturities(dateMaturities)
{
	
	m_nToday	=	m_hUnderlying->GetDateHandle()->GetDate();
	m_nSlice	=	m_dateMaturities.getsize();

	m_t.resize(m_nSlice);
	m_dt.resize(m_nSlice);
	double tprec = 0.0;

	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	
	
	for(int k=0; k<m_nSlice; k++)
	{
		m_t[k]  = hDate->GetYearFraction(m_dateMaturities[k]);	
		m_dt[k] = m_t[k] - tprec;	
		tprec = m_t[k];
	}

	m_params.hParams.resize(m_nSlice);
	m_params.t = m_t ;
}


void heston_term_fitter::calibrateToVarianceSwap()
{
	m_calibVanilla	=	false;
	m_vs.resize(m_nSlice);

	long dateShort = MlEqMaths::Min(m_nToday+90, (m_nToday+m_dateMaturities[0])/2);
	long dateLong  = MlEqMaths::Max(m_nToday+2000, m_dateMaturities[m_nSlice-1] + 365 ) ;

	double varShort = 0., varLong = 0.;
	m_avgKappa = 0.;

	for(int i=0; i<m_nSlice; i++)
	{
		RCPtr<heston_fitter> local_fitter = new heston_fitter(m_hUnderlying, dateShort, m_dateMaturities[i], dateLong);
		m_vs[i] = local_fitter->fComputePrice() / m_t[i] ;

		local_fitter->calibrateToVarianceSwap();
		HestonParameters* local_params = local_fitter->getParams();

		varShort	 += local_params->v0;
		varLong		 += local_params->vas;
		m_avgKappa	 += local_params->kappa;
	}

	varShort /= double(m_nSlice);
	varLong  /= double(m_nSlice);
	m_avgKappa	 /= double(m_nSlice);


	
//////////////////////	
	
	CMatrix NonZeroPartialDerivativesSpecification;
	double FinalTolerance		= 1e-12;/////
	double StoppingTolerance	= 1e-12;////
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;	
	CVector initialGuess(m_nSlice);
	CMatrix xValBounds(m_nSlice,2);
	iVector linearVariableIndex(m_nSlice,0);
	CMatrix ObjectiveBounds;
	CVector x;
	int returnCode;


	for(int i=0;i<m_nSlice;i++)
	{
		initialGuess[i]	= m_avgKappa;
		xValBounds[i][0] = 1e-4;
		xValBounds[i][1] = 10. ;
		m_params.hParams[i].v0	 = varShort;
		m_params.hParams[i].vas	 = varLong;
	}

	initialize(initialGuess,linearVariableIndex,
							xValBounds,ObjectiveBounds,
							InitialTolerance,FinalTolerance,StoppingTolerance,
							NonZeroPartialDerivativesSpecification,outputFlag);

	x = initialGuess ;
	solve(x, 0, returnCode);

	for(int i=0; i<m_nSlice; i++)
	{		
		m_params.hParams[i].kappa = x[i]; //kappa;	x[i] = kappa ;		
	}

	double check;
	ObjectiveFcn(&check, x);
}

void heston_term_fitter::ObjectiveFcnVS(double* gVals,double* xVals)
{
	*gVals = 0.0; 

	double int_kappa = 0.0;
	double varSwap = 0.0;

	double v0 = m_params.hParams[0].v0 ;
	double vas = m_params.hParams[0].vas ;
	double kappa_weight = 0.02 / m_t[m_nSlice-1];	// 1% kappa for 2bp variance swap...

	for(int i=0; i<m_nSlice; i++)
	{
		double	dInt = exp(-int_kappa)*( 1.-exp(-xVals[i]*m_dt[i]) ) / xVals[i];
		varSwap += v0*dInt + vas*(m_dt[i] - dInt);

		double vs_h = varSwap/m_t[i] ;
		vs_h -= m_vs[i] ;
		vs_h = fabs(vs_h);

		double dvs = (vs_h<0.0005)? vs_h*vs_h : vs_h ;

		double dkappa	=	m_avgKappa - xVals[i];
		*gVals += dvs +	kappa_weight * fabs(dkappa) * m_dt[i] ;
		
		int_kappa	  += xVals[i]*m_dt[i];	
	}
}



void heston_term_fitter::calibrate( GVector<MlEqStrikeHandle> vStrikes )
{
	calibrateToVarianceSwap();

//	return;

	double correl = 0.;
	double vovol  = 0.;
	GVector<MlEqStrikeHandle> strikes(2);

//	return;

	for(int slice=0; slice<m_nSlice; slice++)	// calibrate roughly using 0 and -1 norm strikes
	{
		strikes[0] = new MlEqNormalizedStrike(*m_hUnderlying, -1., m_dateMaturities[slice], true, m_hUnderlying->GetDateHandle()) ;
		strikes[1] = new MlEqNormalizedStrike(*m_hUnderlying,  0., m_dateMaturities[slice], true, m_hUnderlying->GetDateHandle()) ;

		calibrateToVanilla( slice, strikes );

		correl += std::min( std::max(m_params.hParams[slice].correl, -1.), 1. ) ;
		vovol  += m_params.hParams[slice].vovol ;
	}

	return;

	vovol  /= double(m_nSlice);		// get average parameters
	correl /= double(m_nSlice);

	for(int slice=0; slice<m_nSlice; slice++)
	{
		m_params.hParams[slice].correl = correl ;
		m_params.hParams[slice].vovol  = vovol ;
	}
	

	for(int slice=0; slice<m_nSlice; slice++)	// now don't calibrate correlation
	{
		calibrateToVanilla( slice, vStrikes, false, correl );	
	}

	return;
}



bool heston_term_fitter::calibrateToVanilla( int slice, GVector<MlEqStrikeHandle> vStrikes, bool calibCorrel, double correl )
{
	m_calibVanilla	=	true;
	m_pricer = new heston_term_pricer(m_params);

	m_pricer->m_nLastSlice	=	m_nCalibSlice	=	slice+1;
	m_calibCorrel = calibCorrel ;
	m_dateMaturity = m_dateMaturities[slice];


	initCalibrationInputs(vStrikes);

	CMatrix NonZeroPartialDerivativesSpecification;

	double FinalTolerance		= 1e-8;
	double StoppingTolerance	= 1e-8;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	int returnCode;

	iVector linearVariableIndex(2,0);
	CMatrix ObjectiveBounds;	
	CVector initialGuess(2);
	CMatrix xValBounds(2,2);
	CVector x;

	if(m_calibCorrel)
	{
		linearVariableIndex.resize(2,0);
		initialGuess.resize(2);
		xValBounds.resize(2,2);

		initialGuess[0]	= -0.2129;//m_cov ;	
		initialGuess[1]	= 0.1954;//m_var ;	
		
		xValBounds[0][0] = -3.;
		xValBounds[0][1] =  3.;
		xValBounds[1][0] =  0.01;
		xValBounds[1][1] =  3.;
	}
	else
	{
		linearVariableIndex.resize(1,0);
		initialGuess.resize(1);
		xValBounds.resize(1,2);

		initialGuess[0]	= 0.25 ;	
	
		xValBounds[0][0] =  0.01;
		xValBounds[0][1] =  3.;

		m_correl = correl;
	}

	initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	x = initialGuess ;	
	solve(x, 0, returnCode);

	double vovol = sqrt( m_pricer->m_vVar[slice] );
	m_params.hParams[slice].vovol  = vovol ;

	if(m_calibCorrel)
	{	
		double correl = m_pricer->m_vCov[slice] / vovol;
		correl = std::min( std::max(correl, -1.), 1. ) ;
		m_params.hParams[slice].correl  = correl ;
	}

	double err = 0.;//heston_term_fitter::checkFcn();

	double check;
	ObjectiveFcn(&check, x);

	return ( err < 0.0010 );	// 10 bp on implied vol	
}





void heston_term_fitter::ObjectiveFcnVanilla(double* gVals,double* xVals)
{
	*gVals = 0.0;

	if(m_calibCorrel)
	{	
		m_pricer->m_vCov[m_nCalibSlice-1] =  xVals[0];
		m_pricer->m_vVar[m_nCalibSlice-1] =  xVals[1];
	}
	else
	{
		m_pricer->m_vVar[m_nCalibSlice-1] =  xVals[0];
		m_pricer->m_vCov[m_nCalibSlice-1] =  sqrt(xVals[0])*m_correl ;
	}

	for(int i=0; i<m_dataSize; i++)
	{
		double price = m_pricer->call_price(m_mat, m_forward, m_strikes[i]);
		price -= m_prices[i] ;
		price /= m_bsVega[i] ;
		*gVals += fabs(price);			
	}
}



sabr_fitter::sabr_fitter( MlEqAssetHandle hUdly, long dateMaturity ):
m_hUnderlying(hUdly),
m_dateMaturity(dateMaturity)
{
}

void sabr_fitter::calibrate( GVector<MlEqStrikeHandle> vStrikes )
{

	initCalibrationInputs(vStrikes);

	CMatrix NonZeroPartialDerivativesSpecification;

	double FinalTolerance		= 1e-8;
	double StoppingTolerance	= 1e-8;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	int returnCode;

	iVector linearVariableIndex(3,0);
	CMatrix ObjectiveBounds;	
	CVector initialGuess(3);
	CMatrix xValBounds(3,2);
	CVector x;

	initialGuess[0]	= m_bsVols[3] ;	
	initialGuess[1]	= 0.5 ;	
	initialGuess[2]	= -0.4 ;	
	
	xValBounds[0][0] =  1e-4;	// v0
	xValBounds[0][1] =  10.;
	xValBounds[1][0] =  0.01;	// vovol
	xValBounds[1][1] =  3.;
	xValBounds[2][0] =  -1.;	// correl
	xValBounds[2][1] =  1.;

	initialize(initialGuess,linearVariableIndex,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);

	x = initialGuess ;	
	solve(x, 0, returnCode);

	m_v0	 = x[0];
	m_vovol	 = x[1];
	m_correl = x[2];


	double check;
	ObjectiveFcn(&check, x);
}



void sabr_fitter::initCalibrationInputs(GVector<MlEqStrikeHandle> vStrikes)
{
	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	
	m_mat	= hDate->GetYearFraction(m_dateMaturity);
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	m_forward =   m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);

	MlEqVolatilityStructureHandle	hVolatilityStructure = m_hUnderlying->GetVolatilityStructure();
	double spot = m_hUnderlying->GetSpot(nToday);	
	
	m_dataSize = vStrikes.getsize();
	if( m_dataSize )
	{
		m_strikes.resize(m_dataSize,0.);
		MlEqStrike stk;	

		for (int i = 0; i < m_dataSize; i++)
		{
			MlEqStrike::convertStrikes(stk,*(vStrikes[i]) );
			m_strikes[i] = stk.m_strike;
		}
	}
	else
	{
		m_dataSize = 6;
		CVector normStrikes(m_dataSize);	
		m_strikes.resize(m_dataSize,0.);

		normStrikes[0] = -2.0;
		normStrikes[1] = -1.0;
		normStrikes[2] = -0.5;
		normStrikes[3] =  0.0;
		normStrikes[4] =  0.5;
		normStrikes[5] =  1.0;

		double vol0 = hVolatilityStructure->getFutureVol(MlEqStrike(m_forward), nToday, m_dateMaturity, spot);

		for(int i=0; i<m_dataSize; i++){
			m_strikes[i] = m_forward* exp(vol0*sqrt(m_mat)* normStrikes[i]) ;	
		}

	}
	
	m_bsVols.resize(m_dataSize,0.);

	for(int i=0; i<m_dataSize; i++)	{
		m_bsVols[i] = hVolatilityStructure->getFutureVol(MlEqStrike(m_strikes[i]), nToday, m_dateMaturity, spot);
	}

}


void sabr_fitter::ObjectiveFcn(double* gVals,double* xVals)
{
	*gVals = 0.0;

	m_v0 = xVals[0];
	m_vovol = xVals[1];
	m_correl = xVals[2];

	for(int i=0; i<m_dataSize; i++)
	{
		double dvol = PatHaagenImpliedVol(m_strikes[i],m_mat,m_forward,m_correl, m_vovol,m_v0);
		dvol -=  m_bsVols[i] ;
		*gVals += fabs(dvol);			
	}
}

void sabr_fitter::getParams(CVector& param)
{
	param.resize(3);

	param[0] = m_v0;
	param[1] = m_vovol;
	param[2] = m_correl;
}




void CStochasticVolMC::initialize( product& deriv, long dateMaturity , int npath,int randomNumberFlag)
{

	int i;

	m_nPaths = npath;
	m_nAssets = 1;

	long nToday = m_hDate->GetDate();
	double T = m_hDate->GetYearFraction(dateMaturity);


	m_nDates = int(T *365./ 1.);	// number of steps
	double dt = T/double(m_nDates);

	m_dt.resize(m_nDates);
	m_t.resize(m_nDates+1);
	m_sqrtdt.resize(m_nDates);

	m_t[0] = 0.0;
	for ( i = 1 ; i <= m_nDates; i++ )
	{
	//	long date = long(nToday + i*dt) ;
		m_t[i] =  i*dt;//m_hDate->GetYearFraction(date);		
	}

	for ( i = 0 ; i < m_nDates; i++ ){
		m_dt[i]		=  m_t[i+1]-m_t[i];	
		m_sqrtdt[i] =  sqrt(m_dt[i]); 
	}

	m_mt = new CMersenneTwister(m_nDates);

	m_numberOfFactors.resize(1);
	m_numberOfFactors[0] = 1;
	m_currentPathValues.resize(m_nAssets,m_nDates+1);




	m_randoms.resize(m_nAssets);
	for ( i = 0 ; i < m_nAssets; i++ ){
		m_randoms[i].resize(m_nDates);
	}


	CMatrix correl(1,1);
	correl[0][0] = 1.;
	m_cholesky.initialize(correl, m_nDates);


// rng stuff

	int dimensionalityFactor = m_nAssets;
	int ndim = dimensionalityFactor*m_nDates;
	int m_randomNumberFlag = randomNumberFlag;

	CUnitcube* D = new CUnitcube(ndim);
			//	select alpha
			//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

	CAlpha* alpha;
	if ( m_randomNumberFlag == 2 ){
				alpha = new CPrimerootalpha;
	}
	else if ( m_randomNumberFlag == 3 ){
				alpha = new CNiederreiteralpha;
	}
	else if ( m_randomNumberFlag == 4 ){
				alpha = new CBakeralpha;
	}
	else if ( m_randomNumberFlag == 5 ){
				alpha = new CPrimerootalpha;
	}
	else{
				throw("randomNumberFlag must be between 2 and 5");
	}

	CArithmeticmeanweights* weights = new CArithmeticmeanweights;// weight
	int d = 1;// d >= r 
	CBakerperiodization* periodization=new CBakerperiodization;
	CNTparameters* par = new CNTparameters(alpha,weights,periodization);
	CNTintegrator* NTintegrator= new CNTintegrator;

	m_randomGenerator = new Diophantine(ndim,m_nPaths,alpha ,
								weights,D,periodization,NTintegrator,par);

}



void CStochasticVolMC::simulate(CMatrix& result,product& deriv)
{

//	generatePaths();

	CMatrix payoffs;
	CMatrix pathArray;
	CVector discounts;

	deriv.setUp(payoffs,*this);

	CMatrix payoffPathArray(m_nAssets, m_nDates);	
	CVector payoffDiscounts(0);

	if ( m_nDates == 0 )
	{
		m_nPaths = 1;
	}

	for (int ipath = 0; ipath < m_nPaths; ipath++ )
	{	
		m_currentPath	= ipath;
		pathArray		= GetPathArray(ipath);

		deriv.payout(payoffs,pathArray,payoffDiscounts,*this);		
		deriv.accruePayout(deriv.m_results,payoffs);
	}
	
	deriv.createResult(deriv.m_results,*this);

	result = deriv.m_results;
}


/*
CMatrix& CSABR_MC::GetPathArray(int ipath)
{
	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);
			
	double prev = 	m_v0*m_v0;
	double average = 0.;

	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		double x = m_vol*m_randoms[0][idate]*m_sqrtdt[idate];
		double next = prev * exp(-m_vol*m_vol*m_dt[idate] + 2.*x) ;		;

		average += 0.5*(prev+next) * m_dt[idate];	
		prev = next;
	}
	
	m_currentPathValues[0][0] = average / m_t[m_nDates-1] ;

	return m_currentPathValues;
}
*/



CMatrix& CSABR_MC::GetPathArray(int ipath)
{	
	double prev = 	m_v0*m_v0;
	double average = 0.;

//	double* g = m_mt->next();
	m_mt->randomGenerator::generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);
	double* g = m_randoms[0].getPtr();

	double* sdt = m_sqrtdt.getPtr();
	double* dt  = m_dt.getPtr();

	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		double ddt = *dt++;
		double x = m_vol * (*g++) * (*sdt++);
		double next = prev * exp(-m_vol*m_vol* ddt + 2.*x) ;		;

		average += 0.5*(prev+next) * ddt ;	
		prev = next;
	}
	
	m_currentPathValues[0][0] = average / m_t[m_nDates-1] ;

	return m_currentPathValues;
}



/*
CMatrix& CHestonMC::GetPathArray(int ipath)
{
	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);			

	double prev = m_v0;
	double average = 0.;

	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		double x = m_randoms[0][idate]*m_sqrtdt[idate];
		double next = prev	+ m_kappa*( m_vas - prev) * m_dt[idate]
							+ m_vol * sqrt(fabs(prev)) * x
							+ 0.25*m_vol*m_vol * (x*x - m_dt[idate]) ;
		

		average += 0.5*(prev+next) * m_dt[idate] ;
		prev = next;
	}	
	
	m_currentPathValues[0][0] = average / m_t[m_nDates-1] ;

	return m_currentPathValues;
}
*/


CMatrix& CHestonMC::GetPathArray(int ipath)
{
//	double* g = m_mt->next();
	m_mt->randomGenerator::generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);
	double* g = m_randoms[0].getPtr();

		
	double* sdt = m_sqrtdt.getPtr();
	double* dt  = m_dt.getPtr();

	double prev = m_v0;
	double average = 0.;

	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		double ddt = *dt++;

		double x = (*g++) * (*sdt++) ;
		double next = prev	+ m_kappa*( m_vas - prev) * ddt
							+ m_vol * sqrt(fabs(prev)) * x
							+ 0.25*m_vol*m_vol * (x*x - ddt) ;
		

		average += 0.5*(prev+next) * ddt ;
		prev = next;
	}	
	
	m_currentPathValues[0][0] = average / m_t[m_nDates-1] ;

	return m_currentPathValues;
}


void CHestonJumpMC::initRandomJump( int rngFlag )
{
	int nMaxJump = 20;	// ???

	m_randomJumps.resize(m_nAssets);
	for (int i = 0 ; i < m_nAssets; i++ ){
		m_randomJumps[i].resize(nMaxJump);
	}

	int dimensionalityFactor = 1;
	int ndim = dimensionalityFactor*nMaxJump;
	int randomNumberFlag = (rngFlag+1)-(rngFlag+1)/3; // so that it is different but between 2 and 5...

	CUnitcube* D = new CUnitcube(ndim);
			//	select alpha
			//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

	CAlpha* alpha;
	if ( randomNumberFlag == 2 ){
				alpha = new CPrimerootalpha;
	}
	else if ( randomNumberFlag == 3 ){
				alpha = new CNiederreiteralpha;
	}
	else if ( randomNumberFlag == 4 ){
				alpha = new CBakeralpha;
	}
	else if ( randomNumberFlag == 5 ){
				alpha = new CPrimerootalpha;
	}
	else{
				throw("randomNumberFlag must be between 2 and 5");
	}

	CArithmeticmeanweights* weights = new CArithmeticmeanweights;// weight
	int d = 1;// d >= r 
	CBakerperiodization* periodization=new CBakerperiodization;
	CNTparameters* par = new CNTparameters(alpha,weights,periodization);
	CNTintegrator* NTintegrator= new CNTintegrator;

	m_randomGeneratorJump = new Diophantine(ndim,m_nPaths,alpha ,
								weights,D,periodization,NTintegrator,par);
}


/*
CMatrix& CHestonJumpMC::GetPathArray(int ipath)
{
	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);			

	double prev = m_v0;
	double average = 0.;

	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		double x = m_randoms[0][idate]*m_sqrtdt[idate];
		double next = prev	+ m_kappa*( m_vas - prev) * m_dt[idate]
							+ m_vol * sqrt(fabs(prev)) * x;
				//			+ 0.25*m_vol*m_vol * (x*x - m_dt[idate]) ;
		

		average += 0.5*(prev+next) * m_dt[idate] ;
		prev = next;
	}
	
	m_randomGeneratorJump->generateUniformRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);			

	double mat = m_t[m_nDates-1];
	double lastJumpTime = 0.;
	int njump = 0;
	double u = m_randoms[0][njump];
	lastJumpTime = -log( u )/m_lambda ;

	while( lastJumpTime < mat && njump < 20)
	{ 
		u = m_randoms[0][++njump] ;
		lastJumpTime += -log(u)/m_lambda ;
	}

	average -= 0.5*( -m_lambda*(exp(m_jump)-1.) + njump*m_jump ) ;
	
	m_currentPathValues[0][0] = average /  mat;

	return m_currentPathValues;
}
*/


CMatrix& CHestonJumpMC::GetPathArray(int ipath)
{
	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);			
	m_randomGeneratorJump->generateUniformRandoms(m_randomJumps,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);	
	
//	double* g = m_mt->next();
	double* g = m_randoms[0].getPtr();
		
	double* sdt = m_sqrtdt.getPtr();
	double* dt  = m_dt.getPtr();



	double prev = m_v0;
	double average = 0.;
	
	double mat = m_t[m_nDates-1];
	double lastJumpTime = 0.;
	int njump = 0;
	double u = m_randomJumps[0][njump];
	lastJumpTime = -log( u )/m_lambda ;


	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		double ddt = *dt++;
		double x = (*g++) * (*sdt++) ;
		double next = prev	+ m_vol * sqrt(fabs(prev)) * x;

		if(lastJumpTime <= m_t[idate] && njump < 20)
		{
			u = m_randomJumps[0][++njump] ;
			lastJumpTime += -log(u)/m_lambda ;

			average -= m_jump * prev * ddt; // average quadratic variation at jump time conditional on V_t (=prev) 
		}

		average += 0.5*(prev+next) * ddt ;
		prev = next;
	}

	average += m_jump*m_jump * double(njump) ;
	
	m_currentPathValues[0][0] = average /  mat;

	return m_currentPathValues;
}

void VarianceCall::init( CVector strikes, double mat, long matDate)
{ 
	m_strikes = strikes; 
	m_mat = mat;m_numberPayouts = strikes.getsize();
	m_payoffDates.resize(1);
	m_payoffDates[0] = matDate;
}

void VarianceCall::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{
	double spot	= pathArray[0][0];

	for(int k=0; k<m_numberPayouts;k++)	
	{
		values[k][0] =	  std::max( (spot - m_strikes[k]), 0.0) ;
	}
}


void VarianceCall::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(m_numberPayouts, 1);
	m_results.resize(m_numberPayouts, 1);
	product::setUp(value,mc);
}







void PayoffTest::init( const CVector& strikes, const CVector& cp, const GVector<long>& matDates)
{ 
	m_strikes = strikes; 
	m_cp = cp; 
	m_numberPayouts = strikes.getsize();

	m_nDates = matDates.getsize()-1;
	m_payoffDates = matDates ;
}

void PayoffTest::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{
	double spot	= pathArray[0][m_nDates];

	for(int k=0; k<m_numberPayouts;k++)	
	{
		values[k][0] =	  std::max( m_cp[k]*(spot - m_strikes[k]), 0.0) ;
	}
}
/*

void PayoffTest::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{
	int nasset = pathArray.rows();
	double spot	= 0.;

	for(int k=0; k<nasset;k++)		{
		double tmp = pathArray[k][0];
				tmp = pathArray[k][1];
		spot += pathArray[k][1] / pathArray[k][0];
	}
	spot /= double(nasset);

	values[0][0] =	  std::max( (spot - 1.), 0.0) ;
}
*/
void PayoffTest::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(m_numberPayouts, 1);
	m_results.resize(m_numberPayouts, 1);
	product::setUp(value,mc);
}





////////////////////////////////////
////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

void straddleHelper::initialize(MlEqInterpolatorHandle initialValue, double spot)
{
	m_spot			= spot;
	m_initialValue	= initialValue;
}

straddleHelper::straddleHelper(MlEqInterpolatorHandle initialValue, double spot)
{
	initialize(  initialValue,  spot);
}

straddleHelper::straddleHelper( double strike,double spot)
{
	initialize( strike, spot );
} 

void straddleHelper::initialize( double strike, double spot )
{
	m_strike	= strike;
	m_spot		= spot;
} 

void straddleHelper::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{
	if( !!m_initialValue )
	{
		for (int ix = 0; ix < nx;ix++ )
		{
			double xx = pde->pde_x(ix);
			u[0][ix] = m_initialValue->getValue( xx );	
		}
	}
	else
	{
		std::vector<double> sx,sv;
		double log_strike = log( m_strike/m_spot );

		for (int ix = 0; ix < nx;ix++ )
		{
			double xx = pde->pde_x(ix);
			u[0][ix] = fabs(m_spot*exp(xx)-m_strike) ;
				
			sx.push_back( xx );
			sv.push_back( u[0][ix] );
			
			if( ix+1 < nx && xx < log_strike - 1e-8 && pde->pde_x(ix+1) > log_strike + 1e-8 )
			{
				sx.push_back( log_strike );
				sv.push_back( 0.0 );
			}
		}

		int n_nx = sx.size();
		CVector x(n_nx), v(n_nx);
		for(int i=0; i<n_nx; ++i)
		{
			x[i] = sx[i];
			v[i] = sv[i];
		}

		int			addTanhWings = 1;	// as in mc model info
		double		cL = 0.01;	
		double		cR = 3.;			
		int			yPower = 2;	

		GVector<CVector> X(1);	X[0] = x;
		GVector<CVector> V(1);	V[0] = v;
		m_initialValue = new MlEqMonotonicSplineInterpolator( X, V, addTanhWings, cL, cR, yPower );
	}

}	
 
void straddleHelper::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;

	lower.a = 1.0;
	lower.b = 0.0;
	lower.c = 0.0;
    lower.d = 0.0;

	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	upper.a = 1.0;
	upper.b = 0.0;
	upper.c = 0.0;
	upper.d = 0.0;	
}



void straddleHelper::finalize( MlEqPdeDriver* pde )
{
	CMatrix u = (pde->pde_full_solution())[0]; // final solution
	int nx = u.cols();

	CVector xv(nx), yv(nx);
	for (int ix = 0; ix < nx; ++ix )
	{
		xv[ix] =  pde->pde_x(ix);
		yv[ix] =  u[0][ix] ;	
	}	

	m_initialValue->reinitialize( xv, yv , 0);	
}











/*

FittedLocalVol::FittedLocalVol(MlEqAssetHandle pAsset)
{
	m_pAsset	= pAsset;
	m_volbump	= 0.;
	m_nToday	= m_pAsset->GetDateHandle()->GetDate();
	m_spot		= pAsset->GetSpot(m_nToday) ;

	m_gridTypeFlag = 0 ;

	m_nFit = 3;
	m_nStrikes.resize(m_nFit);

	m_nStrikes[0] = -1.;
	m_nStrikes[1] =  0.;
	m_nStrikes[2] =  1.;


	m_bsVega.resize(m_nFit);
	m_mktVols.resize(m_nFit);
	m_mktPrices.resize(m_nFit);
	m_cStrikes.resize(m_nFit);

	initInterpolator();

	initialize();

}


void FittedLocalVol::initialize()
{

	int firstMatShift = 92; /////////// #############

	MlEqVolatilityStructureHandle hVS = m_pAsset->GetVolatilityStructure();

	m_nSlice = hVS->getDates().size();

	m_sliceLocalVol.resize(m_nSlice);
	m_sliceDates.resize(m_nSlice);

	m_firstSlice = 0;
	long prevSliceDate = m_nToday;
	double Tprev = 0.0;

	for(int slice=0; slice<m_nSlice; slice++)
	{
		m_sliceLocalVol[slice] = new DupireLocalVol();
		long sliceDate = hVS->getDates()[slice];		

		if( sliceDate >= m_nToday + firstMatShift )
		{
			m_sliceDates[slice] = sliceDate;

			double fwd = m_pAsset->GetForward(sliceDate, false);
			double T = m_pAsset->GetDateHandle()->GetYearFraction(sliceDate);
			double vol0 = m_pAsset->GetVolatility(MlEqStrike(fwd), sliceDate);

			double fnstdev  = (T < 1.)?6.0 - 2.0 * sqrt(T) : 4.0;

			double highSpot = fwd*exp(fnstdev*vol0*sqrt(T));
			double lowSpot  = fwd*exp(- 1.5*fnstdev*vol0*sqrt(T));

			int nt = (T < 1.)?250 - 200 * sqrt(T/1.5) : 50;
			nt *= (T - Tprev);
			int nx = nt;

			GVector<long> timeGrid;
			buildTimeGrid( prevSliceDate, sliceDate, nt, timeGrid);

			m_sliceLocalVol[slice]->initialize( m_pAsset, lowSpot, highSpot, timeGrid, nx, 0.0, 0);

	
			prevSliceDate = sliceDate;
			Tprev = T;
		}
		else
		{
			m_firstSlice++;
			m_sliceDates[slice] = m_nToday;
		}
	}
}
*/
/*
void FittedLocalVol::createTimeGrid( long maturityDate )
{

	int firstMatShift = 92;//// ###############################
	firstMatShift = MlEqMaths::Min( firstMatShift, (maturityDate - m_nToday)/2 );


	MlEqVolatilityStructureHandle hVS = m_pAsset->GetVolatilityStructure();

	m_nSlice = hVS->getDates().size();

	m_sliceLocalVol.resize(m_nSlice);
	m_sliceDates.resize(m_nSlice);

	m_firstSlice = 0;
	long prevSliceDate = m_nToday;
	double Tprev = 0.0;

	std::vector<long> sDates;

	for(int slice=0; slice<m_nSlice; slice++)
	{
		m_sliceLocalVol[slice] = new DupireLocalVol();
		long sliceDate = hVS->getDates()[slice];		

		if( sliceDate >= m_nToday + firstMatShift )
		{			
			m_sliceDates[slice] = MlEqMaths::Min( sliceDate, maturityDate );

			double fwd = m_pAsset->GetForward(sliceDate, false);
			double T = m_pAsset->GetDateHandle()->GetYearFraction(sliceDate);
			double vol0 = m_pAsset->GetVolatility(MlEqStrike(fwd), sliceDate);

			double fnstdev  = (T < 1.)?6.0 - 2.0 * sqrt(T) : 4.0;

			double highSpot = fwd*exp(fnstdev*vol0*sqrt(T));
			double lowSpot  = fwd*exp(- 1.5*fnstdev*vol0*sqrt(T));

			int nt = (T < 1.)?250 - 200 * sqrt(T/1.5) : 50;
			nt *= (T - Tprev);
			int nx = nt;

			GVector<long> timeGrid;
			buildTimeGrid( prevSliceDate, sliceDate, nt, timeGrid);

			m_sliceLocalVol[slice]->initialize( m_pAsset, lowSpot, highSpot, timeGrid, nx, 0.0, 0);


			nt = timeGrid.getsize();
			for(int t=0; t<nt-1; ++t){
				sDates.push_back( timeGrid[t] );
			}

			prevSliceDate = sliceDate;
			Tprev = T;

			if( sliceDate >= maturityDate )
			{
				m_nSlice = MlEqMaths::Min( slice+1, m_nSlice );
				break;
			}
		}
		else
		{
			m_firstSlice++;
			m_sliceDates[slice] = m_nToday;
		}
	}

	sDates.push_back( maturityDate );

	int nt_tot = sDates.size();
	GVector<long> dates(nt_tot);

	for(int t=0; t<nt_tot; ++t){
		dates[t] = sDates[t];
	}

	DupireLocalVol::createTimeGrid( dates );	// now everything is properly "mapped"...
}
*/


void FittedLocalVol::createTimeGrid( const std::vector<long>& payoffDates )
{
	long maturityDate = payoffDates[payoffDates.size()-1];

	int firstMatShift = 92;//// ###############################
	firstMatShift = MlEqMaths::Min( firstMatShift, (maturityDate - m_nToday)/2 );


	MlEqVolatilityStructureHandle hVS = m_pAsset->GetVolatilityStructure();

	m_nSlice = hVS->getDates().size();

	m_sliceLocalVol.resize(m_nSlice);
	m_sliceDates.resize(m_nSlice);

	m_firstSlice = 0;
	long prevSliceDate = m_nToday;
	double Tprev = 0.0;

	std::vector<long> sDates;

	for(int slice=0; slice<m_nSlice; slice++)
	{
		m_sliceLocalVol[slice] = new DupireLocalVol();
		long sliceDate = hVS->getDates()[slice];		

		if( sliceDate >= m_nToday + firstMatShift )
		{			
			m_sliceDates[slice] = MlEqMaths::Min( sliceDate, maturityDate );

			double fwd = m_pAsset->GetForward(sliceDate, false);
			double T = m_pAsset->GetDateHandle()->GetYearFraction(sliceDate);
			double vol0 = m_pAsset->GetVolatility(MlEqStrike(fwd), sliceDate);

			double fnstdev  = (T < 1.)?6.0 - 2.0 * sqrt(T) : 4.0;

			double highSpot = fwd*exp(fnstdev*vol0*sqrt(T));
			double lowSpot  = fwd*exp(- 1.5*fnstdev*vol0*sqrt(T));

			int nsteps = (T < 1.5)?250 - 200 * sqrt(T/1.5) : 50;
			int nx = nsteps * (T - Tprev);

			GVector<long> timeGrid;
			buildTimeGrid( prevSliceDate, sliceDate, nsteps, payoffDates, timeGrid);

			m_sliceLocalVol[slice]->initialize( m_pAsset, lowSpot, highSpot, timeGrid, nx, 0.01, 0);


			int nt = timeGrid.getsize();
			for(int t=0; t<nt-1; ++t){
				sDates.push_back( timeGrid[t] );
			}

			prevSliceDate = sliceDate;
			Tprev = T;

			if( sliceDate >= maturityDate )
			{
				m_nSlice = MlEqMaths::Min( slice+1, m_nSlice );
				break;
			}
		}
		else
		{
			m_firstSlice++;
			m_sliceDates[slice] = m_nToday;
		}
	}

	sDates.push_back( maturityDate );

	int nt_tot = sDates.size();
	GVector<long> dates(nt_tot);

	for(int t=0; t<nt_tot; ++t){
		dates[t] = sDates[t];
	}

	DupireLocalVol::createTimeGrid( dates );	// now everything is properly "mapped"...
}



/*
void FittedLocalVol::initInterpolationGrid( double stdDev, int nx )
{

// build space grid
	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);

	m_yStrikes.resize(m_nTime);
	int sliceMapper = 0;
	
	for(int slice=m_firstSlice; slice<m_nSlice; slice++)
	{
		int nt = m_sliceLocalVol[slice]->m_nTime ;
		if( slice < m_nSlice-1 ){
			nt--;
		}
		GVector<CVector> yStrikes = m_sliceLocalVol[slice]->m_yStrikes ;
		m_nInterpolationSpace = m_sliceLocalVol[slice]->m_nInterpolationSpace ;

		RCPtr<CLocalVolHelper> pLocalVolGrid = m_sliceLocalVol[slice]->m_pLocalVolGrid ;

		for ( int i = 0 ; i < nt; i++ )
		{
			m_yStrikes[ sliceMapper+i ] = yStrikes[i] ;

			m_localVols[sliceMapper+i].resize(m_nInterpolationSpace);
			m_bLocalVols[sliceMapper+i].resize(m_nInterpolationSpace);

			double bump = m_levelBump->getValue( m_times[ sliceMapper+i ] );

			for(int j=0; j<m_nInterpolationSpace; ++j)
			{
				double y = yStrikes[i][j] ;
				
				m_localVols[sliceMapper+i][j] = pLocalVolGrid->getLocalVol( i, y ) + bump ;
				m_bLocalVols[sliceMapper+i][j] = pLocalVolGrid->getBumpedLocalVol( i, y ) + bump;
			}
		}

		sliceMapper += nt;
	}

//	m_forwards.resize( m_nTime );

	m_pLocalVolGrid = new CLocalVolHelper;
	m_pLocalVolGrid->initialize( this );
}

*/


void FittedLocalVol::initInterpolationGrid( double stdDev, int nx )
{
	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);
	m_yStrikes.resize(m_nTime);

	CVector tmp0, tmpK;

	int tt = 0;
	for(int slice=m_firstSlice; slice<m_nSlice; ++slice)
	{
		int nt = m_sliceLocalVol[slice]->m_nTime - 1;
		if( slice == m_nSlice-1 ){	nt++; }

		int nx	= m_sliceLocalVol[slice]->m_nSpace;
		tmp0	= m_sliceLocalVol[slice]->m_fixedStrikes ;
		tmpK	= tmp0;

		for(int t=0; t<nt; ++t)
		{			
			double fwd = m_sliceLocalVol[slice]->m_forwards[t];			
			for(int k=0; k<nx; ++k)	{
				tmpK[k] = log( tmp0[k] / fwd );
			}

			m_yStrikes[tt] = tmpK;
			m_localVols[tt] = m_sliceLocalVol[slice]->getLocalVolGrid()[t] ;
			m_bLocalVols[tt] = m_sliceLocalVol[slice]->getBumpedLocalVolGrid()[t] ;

			tt++;
		}
	}


//	m_bLocalVols = m_localVols;

	m_pLocalVolGrid = new CLocalVolHelper;
	m_pLocalVolGrid->initialize( this );
}

/*
void FittedLocalVol::initialize(MlEqAssetHandle pAsset,
								double lowSpot,
								double highSpot,
								int nx, 
								long maturityDate,
								double vb,
								int gridTypeFlag ) // check the default type
{
	m_pAsset	= pAsset;
	m_nSpace	= nx;
	m_volbump	= vb;
	m_nToday	= m_pAsset->GetDateHandle()->GetDate();
	m_spot		= pAsset->GetSpot(m_nToday) ;

	m_gridTypeFlag = gridTypeFlag ;
//	m_gridTypeFlag = 0 ;

	m_nFit = 3;
	m_nStrikes.resize(m_nFit);

	m_nStrikes[0] = -1.;
	m_nStrikes[1] =  0.;
	m_nStrikes[2] =  1.;

	m_bsVega.resize(m_nFit);
	m_mktVols.resize(m_nFit);
	m_mktPrices.resize(m_nFit);
	m_cStrikes.resize(m_nFit);

	initInterpolator();

	createTimeGrid( maturityDate );

	DupireLocalVol::createSpaceGrid( lowSpot, highSpot, nx );
///	createLocalDrift();	

	fitSurface();

	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);

	if( !m_pLocalVolGrid ){
		initInterpolationGrid();
	}

	DupireLocalVol::createInterpolatedLocalVolGrid();

	m_vbumped = false;
}
*/

void FittedLocalVol::initialize(MlEqAssetHandle pAsset,
								double lowSpot,
								double highSpot,
								int nx, 
								const std::vector<long>& payoffDates,
								double vb,
								int gridTypeFlag ) // check the default type
{
	m_pAsset	= pAsset;
	m_nSpace	= nx;
	m_volbump	= vb;
	m_nToday	= m_pAsset->GetDateHandle()->GetDate();
	m_spot		= pAsset->GetSpot(m_nToday) ;

	m_gridTypeFlag = gridTypeFlag ;


	createTimeGrid( payoffDates );
	

	initMarketData();
	initInterpolator();	
//	fitSurface();

	if( !m_pLocalVolGrid ){
		initInterpolationGrid();
	}

	DupireLocalVol::createSpaceGrid( lowSpot, highSpot, nx );
	DupireLocalVol::createInterpolatedLocalVolGrid();

	m_vbumped = false;

/*
	std::ofstream toto( "C:\\stuff\\localvol.txt");

	for(int t=0; t< m_nTime; ++t)
	{
		toto << m_Dates[t] << '\t';

		for(int k=0; k<m_nSpace; ++k)
		{
			toto << m_localVols[t][k] << '\t';
		}
		toto <<  '\n';
	}
*/	
}

void FittedLocalVol::initMarketData()
{
	
	m_nFit = 3;
	m_nStrikes.resize(m_nFit);

	m_nStrikes[0] = -1.;
	m_nStrikes[1] =  0.;
	m_nStrikes[2] =  1.;

	m_bsVega.resize(m_nFit);
	m_mktVols.resize(m_nFit);
	m_mktPrices.resize(m_nFit);
	m_cStrikes.resize(m_nFit);
}

void FittedLocalVol::initInterpolator()
{
	m_curveBump = new MlEqConstantInterpolator();

	CVector x(1), y(1);
	x[0] = 0.0;
	y[0] = 0.0;
	m_curveBump->initialize(x,y);
}

/*
void FittedLocalVol::buildTimeGrid(long startDate, long endDate, int nt, GVector<long>& timeGrid)
{
	MlEqConstDateHandle	dateh = m_pAsset->GetDateHandle();

	double T = dateh->GetYearFraction(endDate) - dateh->GetYearFraction(startDate) ;
	nt = std::max( nt, 2 );

	double dt = double(endDate-startDate) / (nt-1.);
	long  date = startDate;

	timeGrid.resize(nt);

	for(int i=0; i<nt; ++i)
	{
		date = long( startDate + i*dt );
		timeGrid[i] = date;
	}
	timeGrid[nt-1] = endDate;
}

*/


void FittedLocalVol::buildTimeGrid(long startDate, long endDate, int stepAYear, const std::vector<long>& payoffDates, GVector<long>& timeGrid)
{
	MlEqConstDateHandle	dateh = m_pAsset->GetDateHandle();

	int nsim = payoffDates.size();

	std::vector<long> sDates;
	sDates.push_back(startDate);

	long nextDate = startDate, prevDate = startDate; 
	double Tprev = dateh->GetYearFraction(startDate) ;	
	double T = Tprev;

//  create simulation dates
	for(int i=0; i<nsim; ++i)
	{
		nextDate = MlEqMaths::Min(payoffDates[i], endDate);

		if( nextDate > startDate && nextDate <= endDate )
		{

			T = dateh->GetYearFraction(nextDate) ;			
			int nSteps = int( stepAYear*(T-Tprev) ) + 1;
			double dt = (nextDate - prevDate) / double(nSteps) ;
			if(dt < 1.)	throw "You want too many steps";
			
			for (int j = 1 ; j < nSteps; j++ ) // start at 1 because we already included sim date prev at the previous i 
			{
				long date = long(prevDate + j*dt) ;
				sDates.push_back(date);
			}
			sDates.push_back(nextDate);

			
			prevDate = nextDate;
			Tprev = T;
		}
		if( payoffDates[i] > endDate )	break;
	}

	int nt = sDates.size();
	timeGrid.resize(nt);

	for(int i=0; i<nt; ++i)	{		
		timeGrid[i] = sDates[i];
	}
}





void FittedLocalVol::initialize( int slice )
{
	MlEqVolatilityStructureHandle	hVS = m_pAsset->GetVolatilityStructure();

	double fwd = m_pAsset->GetForward(m_sliceDates[slice], false);
	double T = m_pAsset->GetDateHandle()->GetYearFraction(m_sliceDates[slice]);
	double vol0 = m_pAsset->GetVolatility(MlEqStrike(fwd), m_sliceDates[slice]);
	double df = m_pAsset->GetPayZeroCurve(true)->GetDiscountFactor(m_nToday, m_sliceDates[slice]);

	for(int k=0; k<m_nFit; k++)
	{
		double strike = fwd * exp(vol0*sqrt(T)* m_nStrikes[k]) ;	
		strike = std::min( std::max( strike, 0.8*fwd ), 1.2*fwd );
		double vol  = m_pAsset->GetVolatility(MlEqStrike(strike), m_sliceDates[slice]);
	
		m_mktVols[k]	=	vol ;
		m_mktPrices[k]	=	::Bs(fwd, vol, T, strike, df, 1) + ::Bs(fwd, vol, T, strike, df, -1);	
		m_bsVega[k]		=	(::Bs(fwd, vol+0.01, T, strike, df, 1)+::Bs(fwd, vol+0.01, T, strike, df, -1) - m_mktPrices[k] ) * 100. ;
		m_cStrikes[k]	=	strike ;
	}
	

	CVector x = m_curveBump->getXData();
	CVector y = m_curveBump->getYData();
	int n = x.getsize();
	CVector xt(n+1), yt(n+1);
	for(int i=0; i<n; ++i)
	{
		xt[i] = x[i];
		yt[i] = y[i];
	}
	xt[n] = T;
	yt[n] = 0.0;
	m_curveBump->reinitialize( xt, yt, 0 );

}


void FittedLocalVol::updateInterpolator( double levb )
{

	CVector x = m_curveBump->getXData();
	int n = x.getsize() ;

	if( n < 2 ) throw "Impossible error !";

	CVector y = m_curveBump->getYData();
	y[n-2] = y[n-1] = levb;
	m_curveBump->reinitialize( x, y, 0 );
}



void FittedLocalVol::fitSurface()
{
	for(int slice=m_firstSlice; slice<m_nSlice; slice++)
	{
		initialize( slice );
		FitSlice( slice );
	}	
}


void FittedLocalVol::FitSlice( int slice )
{
	int maxiter = 5;
	double bptol = 3.;//(slice < 10)? 2. + 3.*(slice-m_firstSlice)/double(10-m_firstSlice) : 5. ;
	double tol = bptol * 1e-4;
	double err = 0.0;

	for(int iter=0; iter<maxiter; ++iter)
	{
		err = fitSlice( slice, tol );	
	
		if( fabs(err) < tol )	break;
	}
}


double FittedLocalVol::fitSlice( int slice, double tol )
{
	double dvol = 0.0;

	for(int k=0; k<m_nFit; k++)
	{
		double price = solvePde( slice, m_cStrikes[k] );
		double volbp = (price - m_mktPrices[k]) / m_bsVega[k];
		dvol += volbp;
//		double tmp = m_cStrikes[k] / m_spot;
	}

	dvol /= double(m_nFit);	

	if( fabs(dvol) > tol )
	{
		MlEqConstDateHandle hDate = m_pAsset->GetDateHandle();

		double T =  hDate->GetYearFraction( m_sliceDates[slice] );
		double Tp = hDate->GetYearFraction( m_sliceDates[slice-1] );

		double dvolf = - dvol * sqrt( T/(T - Tp) );

		updateInterpolator( dvolf );
		updateLocalVolSlice( slice );
	}

	return dvol;
}



double FittedLocalVol::solvePde( int slice, double strike )
{
	CVector tmp(1);

	m_straddle = new straddleHelper();
	m_straddle->initialize( strike, m_spot );	

	m_pde = new pdeLocalVol;	

	for(int islice=slice; islice>=m_firstSlice; --islice)
	{
		m_pde->init(m_sliceLocalVol[islice], &(*m_straddle) );
		m_pde->pde_integrate(tmp);
		m_straddle->finalize( &(*m_pde) );
	}

	double price = m_straddle->m_initialValue->getValue( 0.0 );
	return price;
}




void FittedLocalVol::updateLocalVolSlice( int slice )
{
	CVector times = m_sliceLocalVol[slice]->getTimes();
	int nt = times.getsize();
	int nx = m_sliceLocalVol[slice]->getSpotGrid().getsize();

	for(int t=0; t<nt; ++t)
	{
		double bump = m_curveBump->getValue( times[t] );

		for(int k=0; k<nx; ++k)	{
			m_sliceLocalVol[slice]->bumpLocalVol( t, k,  bump );
		}
	}
}











//******************************************* II *******************************

/*
void FittedLocalVol2::initInterpolator()
{
//	m_curveBump = new MlEqCubicSplineInterpolator(); 

	CVector y(m_nFit);

	GVector<CVector> X(1), Y(1);
	X[0] = m_nStrikes;
	Y[0] = y;

	int			addTanhWings = 1;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 2;	
	
	m_curveBump = new MlEqMonotonicSplineInterpolator( 	X, Y, addTanhWings, cL, cR, yPower );

}
*/

/*
void FittedLocalVol2::initInterpolator()
{

	CVector y(m_nFit);
	CVector initG(3,0.0);

	CMatrix xValBounds(3,2);

	for(int k=0; k<3; ++k)
	{
		xValBounds[k][0] = -1.;
		xValBounds[k][1] = 1.;

		initG[k] = 0.1;
	}


	int			addTanhWings = 0;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 1;	

	double final_tol = 1e-4;
	double stop_tol = 1e-8;
	
	m_curveBump = new MlEqAnalyticCurveWithTanhWingEdge(
		m_nStrikes, y, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

}
*/

/*
void FittedLocalVol2::initInterpolator()
{
	CVector xl(m_nFitLeft);
	CVector xr(m_nFitRight);

	CVector yl(m_nFitLeft);
	CVector yr(m_nFitRight);
	CVector initG(3,0.0);

	CMatrix xValBounds(3,2);

	for(int k=0; k<3; ++k)
	{
		xValBounds[k][0] = -1.;
		xValBounds[k][1] = 1.;

		initG[k] = 0.1;
	}


	int			addTanhWings = 0;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 1;	

	double final_tol = 1e-4;
	double stop_tol = 1e-8;
	
	m_leftBump = new MlEqAnalyticCurveWithTanhWingEdge(
		xl, yl, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

	m_rightBump = new MlEqAnalyticCurveWithTanhWingEdge(
		xr, yr, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

}
*/


void FittedLocalVol2::initInterpolator()
{
	CVector xl(m_nFitLeft);
	CVector xr(m_nFitRight);

	for(int k=0; k<m_nFitLeft; ++k){
		xl[k] = k;
	}
	for(int k=0; k<m_nFitRight; ++k){
		xr[k] = k;
	}

	CVector yl(m_nFitLeft);
	CVector yr(m_nFitRight);
	CVector initG(3,0.0);

	CMatrix xValBounds(3,2);

	for(int k=0; k<3; ++k)
	{
		xValBounds[k][0] = -1.;
		xValBounds[k][1] = 1.;

		initG[k] = 0.1;
	}


	int			addTanhWings = 1;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 1;	

	double final_tol = 1e-4;
	double stop_tol = 1e-8;
	
	m_leftBump = new MlEqAnalyticCurveWithTanhWingEdge(
		xl, yl, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

	m_rightBump = new MlEqAnalyticCurveWithTanhWingEdge(
		xr, yr, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

}



/*
void FittedLocalVol2::updateLocalVolSlice( int slice )
{
	int nt = m_sliceLocalVol[slice]->m_nTime;
	int nx = m_sliceLocalVol[slice]->m_nSpace;

	double fwd = m_sliceLocalVol[slice]->m_forwards[nt-1];
	double bump = 0.0;

	for(int k=0; k<nx; ++k)	
	{
		double strike = m_sliceLocalVol[slice]->m_fixedStrikes[k];
		if( strike > fwd ){
			bump = m_rightBump->getValue( strike/fwd );
		}
		else{
			bump = m_leftBump->getValue( strike/fwd );
		}

		for(int t=0; t<nt; ++t)	{
			m_sliceLocalVol[slice]->bumpLocalVol( t, k,  bump );
		}
	}
}
*/

/*
void FittedLocalVol2::updateLocalVolSlice( int slice )
{
	int nt = m_sliceLocalVol[slice]->m_nTime;
	int nx = m_sliceLocalVol[slice]->m_nSpace;

	double fwd = m_sliceLocalVol[slice]->m_forwards[nt-1];
	double bump = 0.0;

	for(int t=0; t<nt; ++t)
	{
		double fwd = m_sliceLocalVol[slice]->m_forwards[t] ;

		for(int k=0; k<nx; ++k)	
		{
			double strike = m_sliceLocalVol[slice]->m_fixedStrikes[k];
			strike = log( strike/fwd );

			if( strike > 0. ){
				bump = m_rightBump->getValue( strike );
			}
			else{
				bump = m_leftBump->getValue( strike );
			}

			m_sliceLocalVol[slice]->bumpLocalVol( t, k,  bump );
		}
	}
}
*/

void FittedLocalVol2::updateLocalVolSlice( int slice )
{
	int nt = m_sliceLocalVol[slice]->m_nTime;
	int nx = m_sliceLocalVol[slice]->m_nSpace;

	double fwd = m_sliceLocalVol[slice]->m_forwards[nt-1];
	double bump = 0.0;

	for(int k=0; k<nx; ++k)	
	{
		double strike = m_sliceLocalVol[slice]->m_fixedStrikes[k];
		if( strike > fwd ){
			bump = m_rightBump->getValue( log(strike/fwd) );
		}
		else{
			bump = m_leftBump->getValue( log(strike/fwd) );
		}

		for(int t=0; t<nt; ++t)	{
			m_sliceLocalVol[slice]->bumpLocalVol( t, k,  bump );
		}
	}
}





double FittedLocalVol2::fitSlice( int slice, double tol )
{
	testFitter(slice) ;
	return 0.;



	double dvol = 0.0;
	double tst = 0.0;
	double avg = 0.0;

	CVector dval( m_nFit, 0.0);

	for(int k=0; k<m_nFit; k++)
	{
		double price = solvePde( slice, m_cStrikes[k] );
		double volbp = (price - m_mktPrices[k]) / m_bsVega[k];
		dvol += fabs(volbp);
		avg += volbp;
		tst = std::max( tst, fabs(volbp) ) ;
		dval[k] = volbp ;
	}

	dvol /= double(m_nFit);	
	avg /= double(m_nFit);	

	if( fabs(tst) > tol )
	{
		MlEqConstDateHandle hDate = m_pAsset->GetDateHandle();

		double T =  hDate->GetYearFraction( m_sliceDates[slice] );
		double Tp = hDate->GetYearFraction( m_sliceDates[slice-1] );
		double wght = sqrt( T/(T - Tp) ) ;

		for(int k=0; k<m_nFit; k++)
		{
			double val = dval[k]; //2.*dval[k] - avg ;
			dval[k] = - wght * val ;
		}

		updateInterpolator( dval );
		updateLocalVolSlice( slice );
	}

	return tst;
	return dvol;
}




void FittedLocalVol2::initialize( int slice )
{
	MlEqVolatilityStructureHandle	hVS = m_pAsset->GetVolatilityStructure();

	double fwd = m_pAsset->GetForward(m_sliceDates[slice], false);
	double T = m_pAsset->GetDateHandle()->GetYearFraction(m_sliceDates[slice]);
	double vol0 = m_pAsset->GetVolatility(MlEqStrike(fwd), m_sliceDates[slice]);
	double df = m_pAsset->GetPayZeroCurve(true)->GetDiscountFactor(m_nToday, m_sliceDates[slice]);

//	m_xdata.resize(m_nFit);
	m_rxdata.resize(m_nFitRight);
	m_lxdata.resize(m_nFitLeft);

	int kl = 0;
	int kr = 0;

	for(int k=0; k<m_nFit; k++)
	{
		double strike = fwd * exp(vol0*sqrt(T)* m_nStrikes[k]) ;	
//		strike = std::min( std::max( strike, 0.8*fwd ), 1.2*fwd );
		double vol  = m_pAsset->GetVolatility(MlEqStrike(strike), m_sliceDates[slice]);
	
		m_mktVols[k]	=	vol ;
		m_mktPrices[k]	=	::Bs(fwd, vol, T, strike, df, 1) + ::Bs(fwd, vol, T, strike, df, -1);	
		m_bsVega[k]		=	(::Bs(fwd, vol+0.01, T, strike, df, 1)+::Bs(fwd, vol+0.01, T, strike, df, -1) - m_mktPrices[k] ) * 100. ;
		m_cStrikes[k]	=	strike ;

//		m_xdata[k]		=	strike / fwd ;
//		double xval		=	strike / fwd ;
//		double xval		=	m_nStrikes[k] ;
		double xval		=	log( strike / fwd );

		if( k < m_nFitLeft ){
			m_lxdata[kl++] = xval ;
		}
		if( k >= m_nFitLeft-1 ){
			m_rxdata[kr++] = xval ;
		}
	}
	

//	CVector bumpVals(m_nFit, 0.0);
//	m_curveBump->reinitialize( m_xdata, bumpVals, 0 );

	CVector bumpVals(m_nFitLeft, 0.0);
	m_leftBump->reinitialize( m_lxdata, bumpVals, 0 );

	bumpVals.resize(m_nFitRight, 0.0);
	m_rightBump->reinitialize( m_rxdata, bumpVals, 0 );


}


void FittedLocalVol2::updateInterpolator( const CVector& dval )
{
//	m_curveBump->reinitialize( m_xdata, dval, 0);;
//	return;

	CVector ldval(m_nFitLeft);
	CVector rdval(m_nFitRight);

	int kl = 0, kr = 0;

	for(int k=0; k<m_nFit; k++)
	{

		if( k < m_nFitLeft ){
			ldval[kl++] = dval[k];
		}
		if( k >= m_nFitLeft-1 ){
			rdval[kr++] = dval[k];
		}
	}

	m_leftBump->reinitialize( m_lxdata, ldval, 0);
	m_rightBump->reinitialize( m_rxdata, rdval, 0);
}



void FittedLocalVol2::initMarketData()
{

	m_nFit = 5;
	m_nStrikes.resize(m_nFit);

	m_nFitLeft = 3;
	m_nFitRight = 3;


	for(int k=0; k<m_nFit; ++k)
	{
		m_nStrikes[k] = -1.5 + 0.6*k ;
	}

	m_bsVega.resize(m_nFit);
	m_mktVols.resize(m_nFit);
	m_mktPrices.resize(m_nFit);
	m_cStrikes.resize(m_nFit);
}


void FittedLocalVol2::testFitter(int slice)
{
	int			addTanhWings = 1;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 1;	

	double final_tol = 1e-4;
	double stop_tol = 1e-8;

	CVector initG(3,0.1);

	CMatrix xValBounds(3,2);

	for(int k=0; k<3; ++k)
	{
		xValBounds[k][0] = -1.;
		xValBounds[k][1] = 1.;
	}
	
	int nt = m_sliceLocalVol[slice]->m_nTime;
	int nx = m_sliceLocalVol[slice]->m_nSpace;

	std::ofstream yoyo( "C:\\stuff\\Intepolator.txt");


	for(int t=0; t<nt; ++t)	
	{
		CVector vol = m_sliceLocalVol[slice]->m_localVols[t] ;

		CVector x = m_sliceLocalVol[slice]->m_fixedStrikes;
		double fwd = m_sliceLocalVol[slice]->m_forwards[t];

		CVector xl(nx/2);
		CVector voll(nx/2);
		for(int k=0; k<nx/2; ++k)
		{
			xl[k] = log( x[k]/fwd );
			voll[k] = vol[k];
		}

		MlEqInterpolatorHandle interp_left 
			= new MlEqAnalyticCurveWithTanhWingEdge(xl, voll, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

		yoyo << m_sliceLocalVol[slice]->m_Dates[t] << '\t';

		for(int k=0; k<nx/2; ++k)
		{
			double dupire = vol[k] ;
			double quadratic = interp_left->getValue( xl[k], 0);

			yoyo << quadratic - dupire << '\t';
		}

		CVector xr(nx-nx/2);
		CVector volr(nx-nx/2);
		for(int k=nx/2; k<nx; ++k)
		{
			xr[k-nx/2] = log( x[k]/fwd );			
			volr[k-nx/2] = vol[k];
		}

		MlEqInterpolatorHandle interp_right 
			= new MlEqAnalyticCurveWithTanhWingEdge(xr, volr, initG, xValBounds, addTanhWings, cL, cR, yPower, final_tol, stop_tol );

		yoyo << m_sliceLocalVol[slice]->m_Dates[t] << '\t';

		for(int k=nx/2; k<nx; ++k)
		{
			double dupire = vol[k] ;
			double quadratic = interp_right->getValue( xr[k-nx/2], 0);

			yoyo << quadratic - dupire << '\t';
		}


		yoyo << '\n';

	}

}