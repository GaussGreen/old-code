
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf hk Model Visitor														 *
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

#include "gpinflation/infhkvisitor.h"

#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal.h"


CC_BEGIN_NAMESPACE( ARM )



/////////////////////////////////////////////////////

///	Class  : ARM_InfHK
///	Routine: ARM_InfHK
///	Returns: void
///	Action : Constructor

/////////////////////////////////////////////////////

ARM_InfHK::ARM_InfHK(	const string & name, 
						const double & dis, 
						const double & dom, 
						const double & eps,
						const double & cen):ARM_InfIntModel(name, dis, dom, eps, cen){

	GaussLegendre_Coefficients c( (int) dis, cen-dom, cen+dom);
	itsWeigt= ARM_GP_VectorPtr ( new ARM_GP_Vector (  (int) dis ) );
	for ( int i = 0; i< (int) dis; i++)
		itsWeigt->Elt(i) =	c.get_weight(i)/sqrt(2.0*PI);
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfHKValue
///	Routine: operator()
///	Returns: double
///	Action : evaluation of the payoff

/////////////////////////////////////////////////////

double ARM_InfHKValue::operator() (			const  ARM_InfPayOffValuePtr	& payOff, 	
											const  double					& lag,
											const  ARM_MAP_Double			& fwd,
											const  ARM_MAP_PairDb			& cor,										
											const  ARM_MAP_VolPar			& vol ){
	ValidatePayOff(payOff);

	if( fwd.size() != vol.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : pb in the NumBiLog model theere incoherenec between inputs");
	
	typedef pair<string,string> Key;

	int dim = fwd.size();

	ARM_MAP_VolPar::const_iterator it;
	ARM_MAP_Double	Par;
	ARM_MAP_Double	Fwd = fwd;

	double			tmp, sum;
	string			stg;


	ARM_GP_VectorPtr	pos = itsModel.GetPosit();
	ARM_GP_VectorPtr	wgt = itsModel.GetWeigt();
	
	int nb = pos->size();
	sum = 0.0;
	if	( dim== 1){
		it = vol.begin();

		ARM_VolParam	vol1 = it->second;
		string			stg1 = it->first;
		double			fwd1 = fwd.find(stg1)->second;
		double			uni1, pos1;

		for ( int i		= 0; i< nb; i++){
			pos1		= pos->Elt(i);		
			uni1		= NormalCDF(pos1);
			Fwd[stg1]	= vol1.itsIdx.CptQua(vol1.itsRes, vol1.itsTen, fwd1,fwd1, tmp);
			sum			+=wgt->Elt(i) * exp( -0.5*pos1*pos1 ) * (*payOff)(lag,Fwd,Par);
		}
	}
	else if ( dim == 2){
		it = vol.begin();

		ARM_VolParam	vol1 = it->second;
		string			stg1 = it->first;
		double			fwd1 = fwd.find(stg1)->second;
		double			pos1, wgt1, nor1, gue1; // nor pour normalisation, gue: guess de la  minimization

		it++;
		ARM_VolParam	vol2 = it->second;
		string			stg2 = it->first;
		double			fwd2 = fwd.find(stg2)->second;
		double			pos2, wgt2, nor2, gue2;

		double			cor12 = cor.find(Key (stg1,stg2))->second;
		double			det   = 1.0-cor12*cor12;

		ARM_GP_Vector	qua1(nb);
		ARM_GP_Vector	qua2(nb);
		ARM_GP_Vector	abs(nb);

		nor1 = 0.0;
		gue1 = fwd1;
		nor2 = 0.0;
		gue2 = fwd2;
FILE* toto = fopen("C:\\momo.txt", "a");

fprintf(toto, "===> PREV SMOOTH \n");

		for ( int i	= 0; i< nb; i++){
			tmp		= NormalCDF(pos->Elt(i));
			abs[i]	= tmp;
			qua1[i]	= vol1.itsIdx.CptQua(vol1.itsRes, vol1.itsTen, fwd1, gue1, tmp);
			qua2[i]	= vol2.itsIdx.CptQua(vol2.itsRes, vol2.itsTen, fwd2, gue2, tmp);
fprintf(toto, "EMU  = %f\n", qua1[i]);
fprintf(toto, "IFRF = %f\n", qua2[i]);
			gue1	= qua1[i];
			gue2	= qua2[i];
		}
fprintf(toto, "===> POST SMOOTH \n");

//		ARM_InfIrIndex::SmoothQua ( abs, qua1 );
//		ARM_InfIrIndex::SmoothQua ( abs, qua2 );

		for (  i = 0; i< nb; i++){		
			tmp		= pos->Elt(i);
			tmp		*=tmp;
			tmp		= exp (-0.5*tmp);
			tmp		*=wgt->Elt(i);

			nor1	+=tmp*qua1[i];
			nor2	+=tmp*qua2[i];
fprintf(toto, "EMU  = %f\n", qua1[i]);
fprintf(toto, "IFRF = %f\n", qua2[i]);
		}
fclose(toto);
		nor1	=	fwd1/nor1;
		nor2	=	fwd2/nor2;
		qua1	*=	nor1;
		qua2	*=	nor2;

		for ( i	= 0; i< nb; i++){
			pos1		= pos->Elt(i);
			wgt1		= wgt->Elt(i);
			Fwd[stg1]	= qua1[i];
			tmp			= pos1*pos1;

			for( int j		= 0; j< nb; j++){
				pos2		= pos->Elt(j);
				wgt2		= wgt->Elt(j);
				Fwd[stg2]	= qua2[j];

				tmp			= pos1*pos1;
				tmp			+=pos2*pos2;
				tmp			-=2.0*cor12*pos1*pos2;
				tmp			/=det;
				tmp			= exp (-0.5*tmp);
				tmp			*=wgt1*wgt2/sqrt(det);

				sum			+= tmp * (*payOff)(lag,Fwd,Par);
			}
		}
	}
	else if ( dim ==3){
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : This code has to be tested");
		it = vol.begin();

		ARM_VolParam	vol1 = it->second;
		string			stg1 = it->first;
		double			fwd1 = fwd.find(stg1)->second;
		double			pos1, wgt1;

		it++;
		ARM_VolParam	vol2 = it->second;
		string			stg2 = it->first;
		double			fwd2 = fwd.find(stg2)->second;
		double			pos2, wgt2;

		it++;
		ARM_VolParam	vol3 = it->second;
		string			stg3 = it->first;
		double			fwd3 = fwd.find(stg3)->second;
		double			pos3, wgt3;

		double cor12= cor.find(Key (stg1,stg2) )->second;
		double cor13= cor.find(Key (stg1,stg3) )->second;
		double cor23= cor.find(Key (stg2,stg3) )->second;
		
		double det = 1.0;

		double c11 = cor12*cor12;
		double c22 = cor13*cor13;
		double c33 = cor23*cor23;

		double c12 = cor13*cor23-cor12;
		double c13 = cor12*cor23-cor13;
		double c23 = cor12*cor13-cor23;

		ARM_GP_Vector	qua1(nb);
		ARM_GP_Vector	qua2(nb);
		ARM_GP_Vector	qua3(nb);
		for ( int i	= 0; i< nb; i++){
			tmp		= NormalCDF(pos->Elt(i)); 
			qua1[i]	= vol1.itsIdx.CptQua(vol1.itsRes, vol1.itsTen, fwd1, fwd1, tmp);
			qua2[i]	= vol2.itsIdx.CptQua(vol2.itsRes, vol2.itsTen, fwd2, fwd2, tmp);
			qua3[i]	= vol3.itsIdx.CptQua(vol3.itsRes, vol3.itsTen, fwd3, fwd3, tmp);
		}


		det -=c11;
		det -=c22;
		det -=c33;
		det +=2.0*cor12*cor13*cor23;

		c11 = 1.0 - c11;
		c22 = 1.0 - c22;
		c33 = 1.0 - c33;
		c12 *=2.0;
		c13 *=2.0;
		c23 *=2.0;
		
		c11 /=det;
		c12 /=det;
		c13 /=det;
		c22 /=det;
		c23 /=det;
		c33 /=det;

		det = sqrt(det);

		for ( i = 0; i< nb; i++){
			pos1		= pos->Elt(i);
			wgt1		= wgt->Elt(i);
			Fwd[stg1]	= qua1[i];

			for( int j		= 0; j< nb; j++){
				pos2		= pos->Elt(j);
				wgt2		= wgt->Elt(j);
				Fwd[stg2]	= qua2[j];

				for ( int k		=0; k< nb ;k++){
					pos3		= pos->Elt(k);
					wgt3		= wgt->Elt(k);
					Fwd[stg3]	= qua3[k];

					tmp			= c11*pos1*pos1;
					tmp			+=c12*pos1*pos2;
					tmp			+=c13*pos1*pos3;
					tmp			+=c22*pos2*pos2;
					tmp			+=c23*pos2*pos3;
					tmp			+=c33*pos3*pos3;

					tmp			= exp(-0.5*tmp);
					tmp			*=wgt1*wgt2*wgt3;
					tmp			/=det;

					sum			+= tmp * (*payOff)(lag,Fwd,Par);
				}
			}
		}
	}
	else {
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : just 3 indexes can be priced");
	}
	return sum;	
}



CC_END_NAMESPACE()





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/



