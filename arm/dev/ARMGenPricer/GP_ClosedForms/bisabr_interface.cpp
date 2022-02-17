/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/bisabr_interface.h"

#include "gpclosedforms/bisabr_digital_spreadoption.h"
#include "gpclosedforms/bisabr_calibration.h"
#include "gpclosedforms/bisabr_spreadoption.h"


CC_BEGIN_NAMESPACE(ARM)


BiSABR_ParameterSet*  Export_BiSABR_CalibrateToSmile(
			ARM_GP_Vector* F1_Vec,ARM_GP_Vector* Alpha1_Vec,ARM_GP_Vector* Beta1_Vec,ARM_GP_Vector* Rho1_Vec,ARM_GP_Vector* Nu1_Vec,
			ARM_GP_Vector* F2_Vec,ARM_GP_Vector* Alpha2_Vec,ARM_GP_Vector* Beta2_Vec,ARM_GP_Vector* Rho2_Vec,ARM_GP_Vector* Nu2_Vec,
			ARM_GP_Vector* K_Vec,ARM_GP_Vector* T_Vec,ARM_GP_Vector* Price_Vec,ARM_GP_Vector* Weight_Vec,
						double rhos_0,double rhov_0,double rhoc12_0,double rhoc21_0,
						double converg_prec, int nbIter_max, double first_step_max,int rhos_flag,int flag)
{
	BiSABR_ParameterSet* setptr;
	if(rhos_flag>0)
	{
		
		setptr= BiSABR_CalibrateToSmile_ConstantRhos(
			F1_Vec,	 Alpha1_Vec,	 Beta1_Vec,	 Rho1_Vec,	 Nu1_Vec,
			F2_Vec,	 Alpha2_Vec,	 Beta2_Vec,	 Rho2_Vec,	 Nu2_Vec,
			K_Vec,	 T_Vec,			 Price_Vec,	 Weight_Vec,
			rhos_0, rhov_0, rhoc12_0, rhoc21_0, 
			converg_prec, nbIter_max, first_step_max,flag);
	}
	else
	{
		setptr= BiSABR_CalibrateToSmile(
			F1_Vec,	 Alpha1_Vec,	 Beta1_Vec,	 Rho1_Vec,	 Nu1_Vec,
			F2_Vec,	 Alpha2_Vec,	 Beta2_Vec,	 Rho2_Vec,	 Nu2_Vec,
			K_Vec,	 T_Vec,			 Price_Vec,  Weight_Vec,
			rhos_0, rhov_0, rhoc12_0, rhoc21_0, converg_prec, nbIter_max, first_step_max,flag);
	}
	
	return setptr;
	
}

void Export_BiSABR_CorrelationEvolution(double rho1,double rho2,double rhos,double rhov,double rhoc12,double rhoc21,
								 double newrho1, double newrho2,
								 double* newrhov,double* newrhoc12,double* newrhoc21)
{
	BiSABR_CorrelationEvolution( rho1, rho2, rhos, rhov, rhoc12, rhoc21,
								  newrho1,  newrho2,
								 newrhov,newrhoc12,newrhoc21);
	return;
}


double Export_BiSABR_Digital_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{


	switch (flag)
			{
			case 0:
				{
						return  BiSABR_Digital_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 1:
				{
						return  BiSABR_Digital_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 2:
				{
						return  BiSABR_Digital_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 3:
				{
						return  BiSABR_Digital_SpreadOption_3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
						break;
				}
			case 4:
				{
					return  BiSABR_Digital_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
						K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}

			case 10:
				{
						return  BiSABR_Digital_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			default:
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_Digital_SpreadOption_PayS1 : flag  bad input :");
						break;
					}
					
			}
}




double Export_BiSABR_Digital_SpreadOption_PayS1(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{
	switch (flag)
			{
			case 0:
				{
					return BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 1:
				{
					return Corrected_BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
						K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 2:
				{
					return Corrected_BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
						K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 4:
				{
					return BiSABR_Digital_SpreadOption_PayS1_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
						K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			case 10:
				{
					return Corrected_BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
						K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
					break;
				}
			default:
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_Digital_SpreadOption_PayS1 : flag  bad input :");
						break;
					}
					
			}
	
}

/// les digitale sont calculées par differentiation de fomule de bas niveau plutot que de haut niveau
/// en faisant appel a Export_BiSABR_SpreadOption(double F1,...,int flag)
double Export_BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{
	switch (flag)
	{
	case 0:
		{
			return BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			break;
		}
	case 1:
		{
			return Corrected_BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			break;
		}
	case 2:
		{
			return Corrected_BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			break;
		}
	case 4:
		{
			return BiSABR_Digital_SpreadOption_PayS2_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			break;
		}
	case 10:
		{
			return Corrected_BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			break;
		}
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_Digital_SpreadOption_PayS2 : flag  bad input :");
			break;
		}
		
	}
	
}


double Export_BiSABR_Digital_SpreadOption_PayS3(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double rhos,double rhov,double rhoc12,double rhoc21,double S3,double sigma3, double rho13, double rho23,double K,double T,int CallPut,int flag)
{
		switch (flag)
		{
		case 0:
			{
				return BiSABR_Digital_SpreadOption_PayS3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					rhos, rhov, rhoc12, rhoc21, S3,sigma3, rho13, rho23,K, T, CallPut);
				break;
			}
		case 1:
			{
				return Corrected_BiSABR_Digital_SpreadOption_PayS3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					rhos, rhov, rhoc12, rhoc21, S3,sigma3, rho13, rho23,K, T, CallPut);
				break;
			}
		case 2:
			{
				return Corrected_BiSABR_Digital_SpreadOption_PayS3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					rhos, rhov, rhoc12, rhoc21, S3,sigma3, rho13, rho23,K, T, CallPut);
				break;
			}
		case 4:
			{
				return BiSABR_Digital_SpreadOption_PayS3_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					rhos, rhov, rhoc12, rhoc21, S3,sigma3, rho13, rho23,K, T, CallPut);
				break;
			}
		case 10:
			{
				return Corrected_BiSABR_Digital_SpreadOption_PayS3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					rhos, rhov, rhoc12, rhoc21, S3,sigma3, rho13, rho23,K, T, CallPut);
				break;
			}
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_Digital_SpreadOption_PayS3 : flag  bad input :");
				break;
			}
			
		}
}




CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/