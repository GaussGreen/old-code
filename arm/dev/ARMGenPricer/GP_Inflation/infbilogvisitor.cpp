
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Bi Log Model Visitor													 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#include "gpinflation/infpayoffvisitor.h"
#include "gpinflation/infhybridpayoffvisitor.h"
#include "gpinflation/infoptionspreadvisitor.h"
#include "gpinflation/infdoubledigitalvisitor.h"

#include "gpinflation/infbilogvisitor.h"

#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/brent.h"
#include "gpnumlib/numfunction.h"

#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal.h"


CC_BEGIN_NAMESPACE( ARM )



/////////////////////////////////////////////////////

///	Class  : ARM_InfNumBiLogValue
///	Routine: CptStrike
///	Returns: double
///	Action : return the dual strike

/////////////////////////////////////////////////////

double	ARM_InfNumBiLogValue::CptStrike (	const string				& str, 
											const double				& lag, 
											const ARM_InfPayOffValuePtr & payOff, 
											const ARM_MAP_Double		& fwd){

	ARM_MAP_Double tmp = fwd;
	double strike  = tmp[str];

	if ( dynamic_cast<ARM_InfHybridPayOffValue*> ( &*payOff)  ) {
		ARM_InfHybridPayOffValue* tmpPayOff = dynamic_cast<ARM_InfHybridPayOffValue*> ( &*payOff);

		tmp[str]=0;
		double res0 = tmpPayOff->CptOpt(lag,tmp);
		tmp[str]=1.0;
		double res1 = tmpPayOff->CptOpt(lag,tmp);
		if ( res1-res0!= 0.0 )	strike	= -res0/(res1-res0);
	}
	return strike;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfNumBiLogValue
///	Routine: CptVol
///	Returns: double
///	Action : return the volatility if computed from call or digital

/////////////////////////////////////////////////////

double call( const double & fwd, const double & str, const double & var ){
	double d1,d2;

	d1 = fwd;
	d1 /=str;
	d1 = log(d1);
	d1 /=var;

	d2 = d1 - 0.5*var;
	d1 +=0.5*var;

	return fwd * NormalCDF(d1) - str * NormalCDF(d2);
}

class CptCallSpread: public ARM_GP::UnaryFunc<double,double>{

public:
	CptCallSpread	(	const double		& fwd, 
						const double		& str,
						const double		& eps,
						const double		& tar): itsFwd(fwd), 
													itsStr(str), 
													itsEps(eps),
													itsTar(tar){ }

	virtual double operator() (double x) const;

private:

	double			itsFwd;
	double			itsStr; 
	double			itsEps;		// epsilon du call spread
	double			itsTar;		// target de la minimization
};



double CptCallSpread::operator() (double x) const{ 
	double tmp;

	tmp = itsFwd;
	tmp /=itsStr;
	tmp = log(tmp);
	tmp /=x;
	tmp -=0.5*x;

	return (call( itsFwd, itsStr, x )-call( itsFwd, itsStr+itsEps, x ))/itsEps-itsTar;		
}


double ARM_InfNumBiLogValue::CptVol	 (	const int						& cri,
										const ARM_InfIrIndex			& idx,
										const double					& res, 
										const double					& ten, 
										const double					& fwd, 
										const double					& str,
										const double					& eps){

	double vol= idx.CptVol( res, ten, fwd, str);
	vol	*= sqrt( fabs(res) );

	if ( cri == DIG) {
		double tar = call( fwd, str, vol );

		vol= idx.CptVol( res, ten, fwd, str);
		vol	*= sqrt( fabs(res) );
		tar -= call( fwd, str+eps, vol );

		tar /= eps;

		CptCallSpread func( fwd, str, eps, tar);
	
		UnaryFuncWithNumDerivative<double> funcWithDeriv(func);
		T_NewtonRaphsonSolverNoThrow<UnaryFuncWithNumDerivative<double> > solver(funcWithDeriv, eps, eps ,eps);

		solver.setInitialGuess(vol);
		vol= solver.Solve();
	}

	return vol;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfNumBiLogValue
///	Routine: operator()
///	Returns: double
///	Action : evaluation of the payoff

/////////////////////////////////////////////////////

double ARM_InfNumBiLogValue::operator() (	const ARM_InfPayOffValuePtr		& payOff, 	
											const double					& lag,
											const ARM_MAP_Double			& fwd,
											const ARM_MAP_PairDb			& cor,							
											const ARM_MAP_VolPar			& vol ){

	ValidatePayOff(payOff);

	if( fwd.size() != vol.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : pb in the NumBiLog model there are incoherences between inputs");
	
	typedef pair<string,string> Key;

	int dim = fwd.size();

	ARM_MAP_VolPar::const_iterator it;

	ARM_MAP_Double	Var;
	ARM_MAP_Double	Par;
	ARM_MAP_Double	Fwd = fwd;

	double			str, tmp, res, ten, sum;
	string			stg;

	ARM_VolParam*	volParam;

	for( it = vol.begin(); it!= vol.end(); it++){
		stg		= it->first;
		str		= CptStrike ( stg, lag, payOff, fwd);
		tmp		= fwd.find(stg)->second;

		volParam= new ARM_VolParam( it->second );

		res		= volParam->itsRes;
		ten		= volParam->itsTen;

		tmp		= CptVol( (*payOff)[stg], volParam->itsIdx, res, ten, tmp, str, itsModel.GetEpsilon() );

		Var.insert(pair<string, double> ( stg, tmp ) );
		if ( volParam ) { delete volParam; volParam=NULL; }
	}


	ARM_GP_VectorPtr	pos = itsModel.GetPosit();
	ARM_GP_VectorPtr	wgt = itsModel.GetWeigt();
	
	int nb = pos->size();
	sum = 0.0;
	if	( dim== 1){
		it = vol.begin();
		string stg1 = it->first;
		double var1 = Var[stg1];
		double fwd1 = fwd.find(stg1)->second * exp( -0.5*var1*var1 );

		for ( int i = 0; i< nb; i++){
			Fwd[stg1] = fwd1*exp( var1*pos->Elt(i) );
			sum += wgt->Elt(i) * (*payOff)(lag,Fwd,Par);
		}
	}
	else if ( dim == 2){
		it = vol.begin();
		double pos1;
		string stg1 = it->first;
		double var1 = Var[stg1];
		double fwd1 = fwd.find(stg1)->second * exp( -0.5*var1*var1);

		it++;
		string stg2 = it->first;
		double var2 = Var[stg2];
		double fwd2 = fwd.find(stg2)->second * exp( -0.5*var2*var2);

		double cor1 = cor.find(Key (stg1,stg2) )->second;
		double cor2 = sqrt(1-cor1*cor1);

		double tmp1;
		for ( int i = 0; i< nb; i++){
			pos1		= pos->Elt(i);
			Fwd[stg1]	= fwd1*exp( var1*pos1 );
			tmp1		= 0.0;

			for( int j	= 0; j< nb; j++){
				Fwd[stg2]= fwd2*exp( var2*(cor1*pos1 + cor2*pos->Elt(j) ) );
				tmp1 += wgt->Elt(j) * (*payOff)(lag,Fwd,Par);
			}
			sum += tmp1 * wgt->Elt(i);
		}
	}
	else if ( dim ==3){
		it = vol.begin();
		double pos1;
		string stg1 = it->first;
		double var1 = Var[stg1];
		double fwd1 = fwd.find(stg1)->second * exp( -0.5*var1*var1 );

		it++;
		double pos2;
		string stg2 = it->first;
		double var2 = Var[stg2];
		double fwd2 = fwd.find(stg2)->second * exp( -0.5*var2*var2 );		
		
		it++;
		string stg3 = it->first;
		double var3 = Var[stg3];
		double fwd3 = fwd.find(stg3)->second * exp( -0.5*var3*var3 );

		double cor21= cor.find(Key (stg1,stg2))->second;
		double cor31= cor.find(Key (stg1,stg3))->second;

		double cor22= sqrt(1-cor21*cor21);
		tmp			= -cor.find(Key(stg1,stg3))->second*cor.find(Key(stg1,stg2))->second; 
		tmp			+=cor.find(Key(stg2,stg3))->second;
		tmp			/=sqrt(1-cor.find(Key (stg1,stg2))->second*cor.find( Key(stg1,stg2) )->second);
		double cor32= tmp;

		tmp			= cor.find(Key (stg1,stg2))->second*cor.find(Key (stg1,stg2))->second;
		tmp			+=cor.find(Key (stg1,stg3))->second*cor.find(Key (stg1,stg3))->second;
		tmp			+=cor.find(Key (stg2,stg3))->second*cor.find(Key (stg2,stg3))->second;
		tmp			= -tmp +1.0+2.0*cor.find(Key (stg1,stg2))->second*cor.find(Key (stg1,stg3))->second*cor.find(Key (stg2,stg3))->second;
		tmp			/=sqrt(1-cor.find(Key (stg1,stg2))->second*cor.find(Key (stg1,stg2))->second);
		double cor33= tmp; 

		double tmp1;
		double tmp2;
		for ( int i = 0; i< nb; i++){
			pos1		= pos->Elt(i);
			Fwd[stg1]	= fwd1*exp(	var1*pos1 );
			tmp2		= 0.0;
			for( int j	= 0; j< nb; j++){
				pos2	= pos->Elt(j);
				Fwd[stg2]= fwd2*exp( var2*(cor21*pos1 + cor22*pos2 )  );
				tmp1 = 0.0;
				for ( int k=0; k< nb ;k++){
					Fwd[stg3]= fwd3*exp( var3*(cor31*pos1 + cor32*pos2 + cor33*pos->Elt(k) ) );
					tmp1 += wgt->Elt(k) * (*payOff)(lag,Fwd,Par);
				}
				tmp2 += tmp1 * wgt->Elt(j);
			}
			sum += tmp2 * wgt->Elt(i);
		}
	}
	else {
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : just 3 indexes can be priced");
	}
	return sum;
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfNumBiLogValue
///	Routine: operator()
///	Returns: double
///	Action : evaluation of the payoff

/////////////////////////////////////////////////////

double ARM_InfBiLogValue::operator() (		const	ARM_InfPayOffValuePtr	& payOff, 	
											const	double					& lag,
											const	ARM_MAP_Double			& fwd,
											const	ARM_MAP_PairDb			& cor,							
											const	ARM_MAP_VolPar			& vol ){

	if( fwd.size() != vol.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : pb in the NumBiLog model there are incoherences between inputs");

	ValidatePayOff(payOff);

	typedef pair<string,string> Key;

	int dim = fwd.size();

	ARM_MAP_VolPar::const_iterator it;

	ARM_MAP_Double	Var;
	ARM_MAP_Double	Fwd = fwd;
	ARM_MAP_Double	Par;

	double			str, tmp, res, ten, sum;
	string			stg;

	ARM_VolParam*	volParam;

	for( it = vol.begin(); it!= vol.end(); it++){
		stg		= it->first;
		str		= CptStrike ( stg, lag, payOff, fwd);
		tmp		= fwd.find(stg)->second;

		volParam= new ARM_VolParam( it->second );

		res		= volParam->itsRes;
		ten		= volParam->itsTen;

		tmp		= CptVol( (*payOff)[stg], volParam->itsIdx, res, ten, tmp, str, itsModel.GetEpsilon());

		Var.insert(pair<string, double> ( stg, tmp ) );
		if ( volParam ) { delete volParam; volParam=NULL; }
	}

	Par.insert(pair<string, double> ( "Ten", ten ) );
	ARM_GP_VectorPtr	pos = itsModel.GetPosit();
	ARM_GP_VectorPtr	wgt = itsModel.GetWeigt();
	
	int nb = pos->size();
	sum = 0.0;
	if	( dim== 1){
		it = vol.begin();
		string stg1 = it->first;
		double var1 = Var[stg1];
		Par.insert(pair<string, double> ( stg, var1 ) );
		sum = (*payOff)(lag,Fwd,Par);
	}
	else if ( dim == 2){
		it = vol.begin();
		double pos1;
		string stg1 = it->first;
		double var1 = Var[stg1];
		double fwd1 = fwd.find(stg1)->second * exp( -0.5*var1*var1);

		it++;
		string stg2 = it->first;
		double var2 = Var[stg2];
		double cor12= cor.find(Key (stg1,stg2) )->second;
		double fwd2 = fwd.find(stg2)->second * exp( -0.5*var2*var2*cor12*cor12);
		Par[stg2] = var2*sqrt(1-cor12*cor12);

		for ( int i = 0; i< nb; i++){
			pos1		= pos->Elt(i);
			Fwd[stg1]	= fwd1*exp( var1*pos1 );
			Fwd[stg2]	= fwd2*exp( var2*cor12*pos1 );

			sum			+=wgt->Elt(i) * (*payOff)(lag,Fwd, Par);
		}
	}
	else {
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : just 2 indexes can be priced");
	}
	return sum;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfBiLogValue
///	Routine: ValidatePayOff
///	Returns: void
///	Action : Validate if the payoff can be priced by this model.

/////////////////////////////////////////////////////

void ARM_InfBiLogValue::ValidatePayOff( const ARM_InfPayOffValuePtr	& payOff){
	if( !dynamic_cast<ARM_InfHybridCapValue*> (&*payOff) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : this model is unable to price other thing than a cap");
}



CC_END_NAMESPACE()




/*--------------------------------------------------------------------------*/
/*---- End of file ----*/



