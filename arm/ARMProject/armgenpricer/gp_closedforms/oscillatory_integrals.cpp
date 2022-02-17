/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file integrals.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/oscillatory_integrals.h"

#include "gpbase/numericconstant.h"


#include "expt.h"
#include <float.h>

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-14
#define ARM_CF_MAXIT 1000

//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////


double OscillatoryIntegral::value()
{
	int j,i;
	double limitinf,limitsup,localsum,sum=0,k;
	GaussLegendre_Coefficients c(legendrePtnb);
	
	switch(runningMode)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case OscillatoryIntegral::MANUAL :
		{
			int n=boundaryList.size();
			
			if (n<2)
			{
				throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value  : bad boundaryList.size");
			}
			if (legendrePtnb<2)
			{
				throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value  : bad legendrePtnb");
			}
			if (Stage_Nb<2)
			{
				throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value  : bad Stage_Nb");
			}

			m_xDetails.resize(0);
			m_intDetails.resize(0);
		
			for(j=0;j<n-1;j++)
			{
				localsum=0.;
				limitinf=boundaryList[j];
				limitsup=boundaryList[j+1];
				double value;
				for(i=0;i<legendrePtnb;i++)
				{
					k=mapper((limitsup-limitinf)/2.*c.get_point(i)+(limitsup+limitinf)/2.);

					value = integrand(k);
					localsum+= value/mapperjacobian(k)*c.get_weight(i);

					m_xDetails.push_back(k);
					m_intDetails.push_back(value);
				}
				sum+=localsum*(limitsup-limitinf)/2.;
				
			}
			break;
		}
	case OscillatoryIntegral::CONTROLED :
		{
			if (legendrePtnb<2)
			{
				throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value  : bad legendrePtnb");
			}
			if (legendrePtnb_FirstStage<2)
			{
				throw  Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value  : bad legendrePtnb for stage 1");
			}
			GaussLegendre_Coefficients c0(legendrePtnb_FirstStage);
	
			int inext;
			if (Stage_Nb<=0)		/// case with automatic determination of the number of stage
			{	
				vector<double> boundaryList;
				boundaryList.push_back(0);
				/// determination of the limit of the first stage
				if(oscillationspeed!=0)
				{
					{
						boundaryList.push_back(inversemapper(2.*ARM_NumericConstants::ARM_PI/oscillationspeed));
					}
				}
				else
				{
					double Default_oscillationspeed=0.01;		// default
					{
						boundaryList.push_back(inversemapper(2.*ARM_NumericConstants::ARM_PI/Default_oscillationspeed));
					}
				}
				/// first stage
				{
					localsum=0.;
					limitinf=boundaryList[0];
					limitsup=boundaryList[1];
					for(i=0;i<legendrePtnb_FirstStage;i++)
					{
						k=mapper((limitsup-limitinf)/2.*c0.get_point(i)+(limitsup+limitinf)/2.);
						
						localsum+= integrand(k)/mapperjacobian(k)*c0.get_weight(i);
					}
					sum+=localsum*(limitsup-limitinf)/2.;
					
				}
				inext=2;
				do
				{
					/// determination of the limit of the next stage
					if(oscillationspeed!=0)
					{
						{
							boundaryList.push_back(inversemapper(2.*inext*ARM_NumericConstants::ARM_PI/oscillationspeed));
						}
					}
					else
					{
						double Default_oscillationspeed=0.01;		// default
						{
							boundaryList.push_back(inversemapper(2.*inext*ARM_NumericConstants::ARM_PI/Default_oscillationspeed));
						}
					}
					/// next stage computation
					{
						localsum=0.;
						limitinf=boundaryList[inext-1];
						limitsup=boundaryList[inext];
						for(i=0;i<legendrePtnb_FirstStage;i++)
						{
							k=mapper((limitsup-limitinf)/2.*c0.get_point(i)+(limitsup+limitinf)/2.);
							
							localsum+= integrand(k)/mapperjacobian(k)*c0.get_weight(i);
						}
						sum+=localsum*(limitsup-limitinf)/2.;
						
					}
					inext++;
				}
				while (fabs(localsum)>SpecifiedPrecision);
			}
			else		// case with number of stage specified
			{
				vector<double> boundaryList(Stage_Nb);
				
				/// Computation of the position of the stages
				if(oscillationspeed!=0)
				{
					for	(i=1;i<Stage_Nb;i++)
					{
						boundaryList[i]=inversemapper(2.*i*ARM_NumericConstants::ARM_PI/oscillationspeed);
					}
				}
				else
				{
					double Default_oscillationspeed=0.01;		// default
					for	(i=1;i<Stage_Nb;i++)
					{
						boundaryList[i]=inversemapper(2.*i*ARM_NumericConstants::ARM_PI/Default_oscillationspeed);
					}
				}
				int n=boundaryList.size();
				
				/// first stage
				{
					localsum=0.;
					limitinf=boundaryList[0];
					limitsup=boundaryList[1];
					for(i=0;i<legendrePtnb_FirstStage;i++)
					{
						k=mapper((limitsup-limitinf)/2.*c0.get_point(i)+(limitsup+limitinf)/2.);
						
						localsum+= integrand(k)/mapperjacobian(k)*c0.get_weight(i);
					}
					sum+=localsum*(limitsup-limitinf)/2.;
					
				}
				
				for(j=1;j<n-1;j++)
				{
					localsum=0.;
					limitinf=boundaryList[j];
					limitsup=boundaryList[j+1];
					for(i=0;i<legendrePtnb;i++)
					{
						k=mapper((limitsup-limitinf)/2.*c.get_point(i)+(limitsup+limitinf)/2.);
						
						localsum+= integrand(k)/mapperjacobian(k)*c.get_weight(i);
					}
					sum+=localsum*(limitsup-limitinf)/2.;
					
				}
			}
			break;
		}
	case OscillatoryIntegral::AUTOADAPTATIF :
		{
			double oscillationdensity,mappingspace_limitsup,mappingspace_limitinf;

			int nbpoint_peroscillation=20;
			int nloop=3;
			mappingspace_limitsup=startpoint;
			limitsup=inversemapper(mappingspace_limitsup);
			for(j=0;j<nloop-1;j++)
			{
				localsum=0.;
				mappingspace_limitinf=mappingspace_limitsup;
				limitinf=limitsup;
				oscillationdensity=fabs(100.*(oscillatorycomponant(mappingspace_limitinf+0.01)-oscillatorycomponant(mappingspace_limitinf))/(2.*ARM_NumericConstants::ARM_PI));
				mappingspace_limitsup=mappingspace_limitinf+(double)legendrePtnb/((double)nbpoint_peroscillation*oscillationdensity);
				limitsup=inversemapper(mappingspace_limitsup);
				for(i=0;i<legendrePtnb;i++)
				{
					k=mapper((limitsup-limitinf)/2.*c.get_point(i)+(limitsup+limitinf)/2.);
					localsum+= integrand(k)/mapperjacobian(k)*c.get_weight(i);
				}
				sum+=localsum*(limitsup-limitinf)/2.;
				
			}


			break;
		}
	
	case NEWADAPTATIF:
		{		
			limitinf	= startpoint;
			limitsup	= getUpperBound();
			
			m_xDetails.resize(0);
			m_intDetails.resize(0);

			double scale = limitsup - limitinf;
			double value;
			for(i = 0; i < legendrePtnb; i++)
			{
				k = limitinf + scale * (0.5 + 0.5 * c.get_point(i));
				value = integrand(k);
				sum += value * 0.5 * c.get_weight(i);

				m_xDetails.push_back(k);
				m_intDetails.push_back(value);
			}
			
			sum *= scale;

		}
		break;

	case NEWMULTADAPTATIF:
		{
			std::vector<double> bounds = getBounds();
			limitinf = startpoint;
			limitsup = bounds[0];

			m_xDetails.resize(0);
			m_intDetails.resize(0);

			for(j = 0; j < (int)bounds.size(); j++)
			{
				

				double scale = limitsup - limitinf;
				double localsum = 0.;
				double value;

				for(i = 0; i < legendrePtnb; i++)
				{
					k = limitinf + scale * (0.5 + 0.5 * c.get_point(i));
					value = integrand(k);
					localsum +=  value * 0.5 * c.get_weight(i);

					m_xDetails.push_back(k);
					m_intDetails.push_back(value);
				}

				localsum *= scale;
				sum += localsum;

				if(j == (int)bounds.size()-1) break;

				limitinf = limitsup;
				limitsup = bounds[j+1];

				
			}
		}
		break;

	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value : incorrect runningMode");
	}
	
	return sum;
	
}


void OscillatoryIntegral::values(std::vector<double>& integrals)
{
	int i,j;
	double limitinf,limitsup,k;
	vector<double> sums(getNbIntegrals(),0.0);
	std::vector<double> integrandVals(getNbIntegrals());
	GaussLegendre_Coefficients c(legendrePtnb);
	
	switch(runningMode)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
	
	case NEWADAPTATIF:
		{		
			limitinf	= startpoint;
			limitsup	= getUpperBound();
			
			double scale = limitsup - limitinf;
			
			for(i = 0; i < legendrePtnb; i++)
			{
				k = limitinf + scale * (0.5 + 0.5 * c.get_point(i));
				integrands(k,integrandVals);
				for(j=0;j<getNbIntegrals();++j)
				{
					sums[j] += integrandVals[j] * 0.5 * c.get_weight(i);
				}
			}
			
			for(j=0;j<getNbIntegrals();++j)
				integrals[j] = sums[j]*scale;
		}
		break;

	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OscillatoryIntegral::value : incorrect runningMode");
	}
}

double OscillatoryIntegral::getUpperBound()
{
	double point = startSearch;
	double prevpoint = point;
	double step = startSearch, eps = 1e-8;

	double fp = integrand(point);
	double prevfp = fp;

	while(fabs(fp) > eps || fabs(prevfp) > eps)
	{
		prevfp = fp;
		prevpoint = point;
		point += step;
		fp = integrand(point);
		if(_isnan(fp))
		{
			point	-= step;
			step	/= 2.;
			fp		= prevfp;
		}
		
		if(fabs(fp) < eps && fabs(prevfp) < eps)
		{
			point = 0.5*(prevpoint + point);
			fp = integrand(point);
		}

		if(point >= 150.) break;
	}

	return point;
}

std::vector<double> OscillatoryIntegral::getBounds()
{
	std::vector<double> bounds;
	double point = startSearch;
	double prevpoint = point;
	double step = startSearch, eps = 1e-8;

	double fp = integrand(point);
	double prevfp = fp;

	while(fabs(fp) > eps || fabs(prevfp) > eps)
	{
		if(fp < 0. && (int)bounds.size() < Stage_Nb) bounds.push_back(prevpoint);
			
		prevfp = fp;
		prevpoint = point;
		point += step;
		fp = integrand(point);
		if(_isnan(fp))
		{
			point	-= step;
			step	/= 2.;
			fp		= prevfp;
		}
		
		if(fabs(fp) < eps && fabs(prevfp) < eps)
		{
			point = 0.5*(prevpoint + point);
			fp = integrand(point);
		}

		if(point >= 150.) break;
	}

	bounds.push_back(point);
	return bounds;
}


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
