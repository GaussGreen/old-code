
#include "Magnet/Magnet.h"
#include "General/General.h"

#include "Convolution.h"
#include "Quadrature.h"
#include "imsl.h"
#include "imsls.h"
#include "Util.h"


using namespace CM;


// Finds the unconditional default probability in Andersen and Sidenius (2003) 
// random factor model given mean, variance, high and low factors as well as threshold
// cf. p22
double UncondSurvProb(
		const double   c, // default threshold
		const double   mean,
		const double   stddev,
		const double   betaH,
		const double   betaL,
		const double   threshold) // market factor threshold
{ 
	double prob;

	double fH=sqrt(stddev*stddev+betaH*betaH);
	double fL=sqrt(stddev*stddev+betaL*betaL);

	prob=imsls_d_bivariate_normal_cdf((c-mean)/fH,threshold,betaH/fH)+GtoNormalCum((c-mean)/fL);
	prob=prob-imsls_d_bivariate_normal_cdf((c-mean)/fL,threshold,betaL/fL);

	return prob;
}


// Where is this used? 
double UncondSurvProbGen(
		double   c,
		double   mean,
		double   stddev,
		double   beta,
		double   DispGen,
		double   Disp1,
		double   Disp2,
		double   theta1,
		double   theta2)
{ 
    // put theta and disp in order
	if (theta1>theta2)
	{ 
		double temp=theta1;
	    theta1=theta2;
	    theta2=temp;
		double temp2=Disp1;
		Disp1=Disp2;
		Disp2=temp2;
	}	
	
	double s1=beta+DispGen+Disp1+Disp2;
	double s2=beta+DispGen-Disp1+Disp2;
	double s3=beta+DispGen-Disp1-Disp2;
	double s4=beta-DispGen+Disp1+Disp2;
	double s5=beta-DispGen-Disp1+Disp2;
	double s6=beta-DispGen-Disp1-Disp2;
	double f1=sqrt(stddev*stddev+s1*s1);
	double f2=sqrt(stddev*stddev+s2*s2);
	double f3=sqrt(stddev*stddev+s3*s3);
	double f4=sqrt(stddev*stddev+s4*s4);
	double f5=sqrt(stddev*stddev+s5*s5);
	double f6=sqrt(stddev*stddev+s6*s6);


	double prob=imsls_d_bivariate_normal_cdf((c-mean)/f1,theta1,s1/f1)+imsls_d_bivariate_normal_cdf((c-mean)/f2,theta2,s2/f2)-imsls_d_bivariate_normal_cdf((c-mean)/f2,theta1,s2/f2);
	prob=prob+GtoNormalCum((c-mean)/f3)-imsls_d_bivariate_normal_cdf((c-mean)/f3,theta2,s3/f3);
	
	prob=0.5*prob+0.5*imsls_d_bivariate_normal_cdf((c-mean)/f4,theta1,s4/f4)+0.5*imsls_d_bivariate_normal_cdf((c-mean)/f5,theta2,s5/f5)-0.5*imsls_d_bivariate_normal_cdf((c-mean)/f5,theta1,s5/f5);
	prob=prob+0.5*(GtoNormalCum((c-mean)/f6)-imsls_d_bivariate_normal_cdf((c-mean)/f6,theta2,s6/f6));

	return prob;
}


// Calc unconditional default prob in poly regime
double UncondSurvProbPolExp(
		double   c,
		double   mean,
		double   stddev,
		double   a0,
		double   a1,
		double   a2,
		double   aexp,
		double   kexp)
{ 
	double weights[40];
	double points[40];

	imsl_d_gauss_quad_rule(40,weights, points, IMSL_HERMITE,0);

	double prob=0;

	double z;
	for (int j=0; j<40; j++)
	{	 
	  z=points[j]*sqrt(2);
	  prob=prob+GtoNormalCum((c-mean-a0*z-a1*z*z-a2*z*z*z-aexp*exp(-kexp*z)*z)/stddev)*weights[j];
	}

	prob=prob/sqrt(3.141592653);
	
	return prob;
}



// Calc unconditional default prob in discrete regime
double UncondSurvProbDiscrete(
		double			c,
		double		    mean,
		double			stddev,
		Array<double>	&prob,
		Array<double>	&theta)

{	
	int NoPoints=theta.size();

	double uncondProb=0;
    for (int s=0; s<(NoPoints); s++)
	{
        // prob[s] is the prob of market factor = s
		uncondProb=uncondProb+GtoNormalCum((c-theta[s]-mean)/stddev)*prob[s];
	}
		
	return uncondProb;
}




// Calc unconditional default prob in general n factor regime
double UncondSurvProbNfact(
		double			c,
		double			mean,
		double			stddev,
		Array<double>   beta,
		Array<double>	theta)
{ 
    int NoBeta=beta.size();
	
	double fi=sqrt(stddev*stddev+beta[0]*beta[0]);
	double prob=imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[0],beta[0]/fi);	
	 
    for (int s=1; s<(NoBeta-1); s++)
	{
		fi=sqrt(stddev*stddev+beta[s]*beta[s]);
		prob=prob+imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[s],beta[s]/fi)
                    -imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[s-1],beta[s]/fi);
	}
	
	fi=sqrt(stddev*stddev+beta[NoBeta-1]*beta[NoBeta-1]);
	prob=prob+GtoNormalCum((c-mean)/fi)
                -imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[NoBeta-2],beta[NoBeta-1]/fi);
			 
	return prob;
}


// Calc unconditional default prob in general n factor regime
// assume market factor follows t distribution, to do
double UncondSurvProbNfactT(
		double			c,
		double			mean,
		double			stddev,
		Array<double>   beta,
		Array<double>	theta,
        int             dof)
{ 
    int NoBeta=beta.size();
	
	double fi=sqrt(stddev*stddev+beta[0]*beta[0]);
	double prob=imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[0],beta[0]/fi);	
	 
    for (int s=1; s<(NoBeta-1); s++)
	{
		fi=sqrt(stddev*stddev+beta[s]*beta[s]);
		prob=prob+imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[s],beta[s]/fi)
                    -imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[s-1],beta[s]/fi);
	}
	
	fi=sqrt(stddev*stddev+beta[NoBeta-1]*beta[NoBeta-1]);
	prob=prob+GtoNormalCum((c-mean)/fi)
                -imsls_d_bivariate_normal_cdf((c-mean)/fi,theta[NoBeta-2],beta[NoBeta-1]/fi);
			 
	return prob;
}


// Given an array of betas and recoveries with their corresponding thresholds
// combines it into one array with piecewise constant betas and recoveries
// (recoveries and beta thresholds are combined and sorted)

Matrix<double> OrderBetaRecTh(
		Array<double>   &Beta, // N betas
		Array<double>	&theta, // N-1 thetas
		Array<double>   &Recovery, // likely 3 
		Array<double>	&thetaR) // likely 2 
{ 
	int NoBeta=Beta.size();
	int NoRec=Recovery.size();

	Array<double> BetaAdj(NoBeta+NoRec-1);
	Array<double> RecAdj(NoBeta+NoRec-1);
	Array<double> ThetaAdj(NoBeta+NoRec-2);

	for (int i=0; i<NoBeta; i++)
	{
		BetaAdj[i]=Beta[i];
		RecAdj[i]=Recovery[0];
		if (i<(NoBeta-1))
		{
			ThetaAdj[i]=theta[i];
		}
	}

	for (int h=0; h<(NoRec-1); h++)
	{
        int index=0;
        while ((ThetaAdj[index]<thetaR[h])&&(index<(NoBeta-1+h)))
		  index=index+1;
		
		for(int l=NoBeta-1+h; l>=index; l--)
		{
            BetaAdj[l+1]=BetaAdj[l];
		  
            double Ah=BetaAdj[l+1];
            if (l<(NoBeta-1+h))
            { 
                ThetaAdj[l+1]=ThetaAdj[l];
            }				
		
            RecAdj[l+1]=Recovery[h+1];
		}

        BetaAdj[index]=BetaAdj[index+1]; 	
        ThetaAdj[index]=thetaR[h];	
	}	

    Matrix<double> param(NoBeta+NoRec-1,3);
    for (int k=0; k<NoBeta+NoRec-1; k++)
    {
        param[k][0]=BetaAdj[k];
        param[k][1]=RecAdj[k];
		
        if (k<NoBeta+NoRec-2)
        {
            param[k][2]=ThetaAdj[k];
        }
    }

    dump(param, -1, -1, "sorted");

    return param;
}	


// Calc average recovery from input params in N random factor, N random recovery regime
double UncondSurvProbNBetaNRec(
		double			c,
		double			mean,
		double			stddev,
		Matrix<double>  &Param)
{	
	int NoBeta=Param.rowSize();
			
    // param[0][0] = beta[0], param[0][1] = rec[0], param[0][2] = theta[0]
	double fi=sqrt(stddev*stddev+Param[0][0]*Param[0][0]);
    double prob=imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[0][2],Param[0][0]/fi);
	double rec=imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[0][2],Param[0][0]/fi)*Param[0][1];
	 
    for (int s=1; s<(NoBeta-1); s++)
	{
		fi=sqrt(stddev*stddev+Param[s][0]*Param[s][0]);
		prob=prob+imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[s][2],Param[s][0]/fi)
                    -imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[s-1][2],Param[s][0]/fi);
		rec=rec+Param[s][1]*(imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[s][2],Param[s][0]/fi)
                                -imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[s-1][2],Param[s][0]/fi));
	}
	
	fi=sqrt(stddev*stddev+Param[NoBeta-1][0]*Param[NoBeta-1][0]);
	prob=prob+GtoNormalCum((c-mean)/fi)
            -imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[NoBeta-2][2],Param[NoBeta-1][0]/fi);
	rec=rec+Param[NoBeta-1][1]*(GtoNormalCum((c-mean)/fi)
                    -imsls_d_bivariate_normal_cdf((c-mean)/fi,Param[NoBeta-2][2],Param[NoBeta-1][0]/fi));
			 
	double avgRecgivenDefault=rec/prob;

    dump(rec, "rec");
    dump(prob, "prob");
    dump(avgRecgivenDefault, "avgRecgivenDefault");

	return avgRecgivenDefault;
}


// Calculates loss dist for a portfolio with homogeneous correlation
// and default probabilities
// Integration over the market factor is done by Gaussian Quadrature

Array<double> LossDistConstBetaHomogeneousSpread(
    const int       NoNames,
    double			Beta,
    double			DefaultProb)
{ 
	
	int MaxNoUnits=NoNames;
	
	Array<double> Loss(MaxNoUnits+1);
	Array<double> LossProb(MaxNoUnits+1);
	  
	double ci;
    GtoNormalCumInv(DefaultProb,&ci); // get default threshold
	
	double weights[40];
	double points[40];

    // finds quadrature weights and points using IMSL routine
    // points essentially stores market factor
	imsl_d_gauss_quad_rule(40,weights, points, IMSL_HERMITE,0);

	Matrix<double> CondLoss(MaxNoUnits+1, 40);
	Array<double> LossTemp(MaxNoUnits+1);

	for (int i=0; i<40; i++)
	{		
		double condDefaultProbi;
	
        // default prob of individual name given current market condition
		condDefaultProbi=GtoNormalCum((ci-Beta*points[i]*sqrt(2))/(sqrt(1-Beta*Beta)));	

        // since all the names have the same default probabilities no Convolution procedure
        // is needed to calculate conditional loss distribution
		 		
		for (int p=0; p<=MaxNoUnits; p++)
		{			
            // gamma(N) = factorial(N-1)
			double comb=imsl_d_gamma(MaxNoUnits+1)/imsl_d_gamma(p+1)/imsl_d_gamma(MaxNoUnits-p+1);
            // prob of having p losses
			LossTemp[p]= comb*pow(condDefaultProbi,p)*pow(1-condDefaultProbi,MaxNoUnits-p);
		}

		for (int k=0; k<=MaxNoUnits; k++)
		{ 
			CondLoss[k][i] = LossTemp[k];
		}		
	}

	for (int j=0; j<=MaxNoUnits; j++)
	{	
		Loss[j] = 0;
		for (int i=0; i<40; i++)
		{ 
			Loss[j] = Loss[j]+weights[i]*CondLoss[j][i]/sqrt(3.141592653);	
		}
	}

	return Loss;
}





// Calculates loss dist for a portfolio with homogeneous correlation
// Loss given default, and default probabilities can differ among
// names
// A standard convolution routine is implemented and numerical integration
// is done by Gaussian Quadrature

Array<double> LossDistConstBeta(
	 const int      NoNames,
	 double			Beta,
	 Array<int>     &LossGivenDefault,
	 Array<double>  &DefaultProb)
{ 
	// finds the loss unit taken as the Greatest Common Divisor over all names of losses given default
	// and the total Maximum Number of Loss Units as the maximum total loss that could be suffered
	// measured in loss units
	
	Array<int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);

	int MaxNoUnits=LossParam[1];	
	
	Array<double> Loss(MaxNoUnits+1);
	Array<double> ci(NoNames);
	Array<double> LossProb(MaxNoUnits+1);

    // calc default threshold for each name
	for (int k=0; k<NoNames; k++)
        GtoNormalCumInv(DefaultProb[k],&ci[k]);

	double weights[40];
	double points[40];

	imsl_d_gauss_quad_rule(40,weights, points, IMSL_HERMITE,0);

	Matrix<double> CondLoss(MaxNoUnits+1, 40);
	Array<double> LossTemp(MaxNoUnits+1);

	for (int i=0; i<40; i++)
	{
		Array<double> condDefaultProbi(NoNames);

		for (int j=0; j<NoNames; j++)	
            condDefaultProbi[j]=GtoNormalCum((ci[j]-Beta*points[i]*sqrt(2))/(sqrt(1-Beta*Beta)));	
	
		// calculates conditional distribution by a standard convolution 
		// procedure
		LossTemp = CondLossDist(NoNames, LossGivenDefault, condDefaultProbi);

		for (int k=0; k<=MaxNoUnits; k++)
		{ 
			CondLoss[k][i] = LossTemp[k];
		}
	}

	for (int j=0; j<=MaxNoUnits; j++)
	{		
		Loss[j] = 0;
		for (int i=0; i<40; i++)
		{ 
			Loss[j] = Loss[j]+weights[i]*CondLoss[j][i]/sqrt(3.141592653);
		}
	}

	return Loss;
}



// Same as previous code with beta allowed to differ
// between names

Array<double> LossDist(
	 const int      NoNames,
	 Array<double>  &Correlation,
	 Array<int>     &LossGivenDefault,
	 Array<double>  &DefaultProb)
{ 
	Array<int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
	int MaxNoUnits=LossParam[1];

	Array<double> Loss(MaxNoUnits+1);	  	
	Array<double> LossProb(MaxNoUnits+1);

	Array<double> ci(NoNames);
	for (int k=0; k<NoNames; k++)
        GtoNormalCumInv(DefaultProb[k],&ci[k]);

	double weights[40];
	double points[40];

	imsl_d_gauss_quad_rule(40,weights, points, IMSL_HERMITE,0);


	Matrix<double> CondLoss(MaxNoUnits+1, 40);
	Array<double> LossTemp(MaxNoUnits+1);

	for (int i=0; i<40; i++)
	{
		Array<double> condDefaultProbi(NoNames);

		for (int j=0; j<NoNames; j++)
		{
            condDefaultProbi[j]=GtoNormalCum((ci[j]-Correlation[j]*points[i]*sqrt(2))/(sqrt(1-Correlation[j]*Correlation[j])));
		}
		 
		LossTemp = CondLossDist(NoNames, LossGivenDefault, condDefaultProbi);

		for (int k=0; k<=MaxNoUnits; k++)
		{ 
			CondLoss[k][i] = LossTemp[k];
		}		
	}


	for (int j=0; j<=MaxNoUnits; j++)
	{
		Loss[j] = 0;
		for (int i=0; i<40; i++)
		{ 
			Loss[j] = Loss[j]+weights[i]*CondLoss[j][i]/sqrt(3.141592653);
		}
	}

	return Loss;
}



// Calc loss distribution in poly regime, homogenous with different default prob

Array<double> LossDistRandomFactorPolExp(
	 const int	   NoNames,
	 double		   a0,
	 double		   a1,
	 double		   a2,
	 double		   aexp,
	 double		   kexp,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb)
{ 
	
	Array<int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
	int MaxNoUnits=LossParam[1];
	
	Array<double> Loss(MaxNoUnits+1);	  	
	Array<double> LossProb(MaxNoUnits+1);

	// calc parameters		 
	double UncondMean = -a1+kexp*aexp*exp(kexp*kexp/2); 
	double V=a0*a0+3*a1*a1+15*a2*a2+aexp*aexp*exp(2*kexp*kexp)*(1+4*kexp*kexp)+6*a0*a2+2*a0*aexp*exp(kexp*kexp/2)*(1+kexp*kexp);
	V=V-2*a1*aexp*exp(kexp*kexp/2)*(3*kexp+kexp*kexp*kexp)+2*a2*aexp*exp(kexp*kexp/2)*(kexp*kexp*kexp*kexp+6*kexp*kexp+3);
	V=V-UncondMean*UncondMean;

	if (V>1)		 
    {
        for (int j=0; j<MaxNoUnits+1; j++)
            Loss[j]=-1;
        return Loss;
    }
    double UncondStdDev=sqrt(1-V);

    // calculate default threshold to match given default prob
    double ci;
    GtoNormalCumInv(DefaultProb[0],&ci);
    {
        double clow;
        double chigh;
        double cmid; 
	    const double alphafact=4;
	    
	    if (ci>0)
        {
            clow=ci/alphafact;
            chigh=ci*alphafact;
            cmid=(clow+chigh)/2;
        } else {
            clow=ci*alphafact;
            chigh=ci/alphafact;
            cmid=(clow+chigh)/2;
        }

        // how far are we away from result
        double zero=UncondSurvProbPolExp(cmid,UncondMean,UncondStdDev,a0,a1,a2,aexp,kexp)-DefaultProb[0]; 
        while (fabs(zero)>0.0000001)
		{
            if (zero>0) // c too big
            {	
			    chigh=cmid;
			    cmid=(chigh+clow)/2;
            } else {
			    clow=cmid;
			    cmid=(chigh+clow)/2;
            }
            zero=UncondSurvProbPolExp(cmid,UncondMean,UncondStdDev,a0,a1,a2,aexp,kexp)-DefaultProb[0];
		}
        ci=cmid; // found c
    }
	
	Matrix<double> CondLoss(MaxNoUnits+1, 40);
	Array<double> LossTemp(MaxNoUnits+1);

    double weights[40];
	double points[40];
	imsl_d_gauss_quad_rule(40,weights, points, IMSL_HERMITE,0);
	
	for (int i=0; i<40; i++)
	{		
		Array<double> condDefaultProbi(NoNames);

        double z;
		for (int j=0; j<NoNames; j++)
		{	
		    z=points[i]*sqrt(2);
  		    condDefaultProbi[j]=GtoNormalCum((ci-(a0+a1*z+a2*z*z+aexp*exp(-kexp*z))*z-UncondMean)/UncondStdDev);	
		}
		 
		LossTemp = CondLossDist(NoNames, LossGivenDefault, condDefaultProbi);
		
		for (int k=0; k<=MaxNoUnits; k++)
		{ 
			CondLoss[k][i] = LossTemp[k];
		}	
	}

	for (int j=0; j<=MaxNoUnits; j++)
	{
		Loss[j] = 0;
		for (int i=0; i<40; i++)
		{ 
			Loss[j] = Loss[j]+weights[i]*CondLoss[j][i]/sqrt(3.141592653);
		}
	}

	return Loss;
}


// Calc distribution in discrete regime using plain Simpson Integration

Array<double> LossDistRandomDiscrete(
	 const int	   NoNames,
	 Array<double> probMarket,
	 Array<double> theta,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb)
{ 
	
	Array<int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
	
	int NoPoints=theta.size();	
	int MaxNoUnits=LossParam[1];
	
	Array<double> Loss(MaxNoUnits+1);	  
	Array<double> ci(NoNames);
	
	Array<double> LossProb(MaxNoUnits+1);

	
	double z;

		

	for (int k=0; k<NoNames; k++)

	{ 
      GtoNormalCumInv(DefaultProb[k],&z);
	  ci[k]=z;
	}
 

	double UncondMean=0;

	double V=0; 

	for (int j=0; j<NoPoints; j++)

	{

		UncondMean=UncondMean+probMarket[j]*theta[j];
		V=V+probMarket[j]*theta[j]*theta[j];
	}


	 UncondMean=-UncondMean;
	 V=V-UncondMean*UncondMean;
	 
 	 if (V>1)
		 
	 {
			 for (int j=0; j<MaxNoUnits+1; j++)
		 {
			 Loss[j]=-1;
		 }
	 return Loss;
	 }
	 

	 double UncondStdDev=sqrt(1-V);


	Array<double> clow(NoNames);
	Array<double> chigh(NoNames);
	Array<double> cmid(NoNames);
	 double alpha=1.5;
	 double zero;
	

for (k=0; k<NoNames; k++)

	{ 
	  if (ci[k]>0)
		
	  {
		clow[k]=ci[k]/alpha;
		chigh[k]=ci[k]*alpha;
		cmid[k]=(clow[k]+chigh[k])/2;

		double zeroH=UncondSurvProbDiscrete(chigh[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k]; 		  

        	while (zeroH<0)

				{
    				alpha=1.2*alpha;
					chigh[k]=chigh[k]*alpha;
					zeroH=UncondSurvProbDiscrete(chigh[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
		  		}
			
    	double zeroL=UncondSurvProbDiscrete(clow[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
	
    	  alpha=1.5; 

				  
			if (zeroL>0)

			{
    	    	alpha=1.2*alpha;
				clow[k]=clow[k]/alpha;
				zeroL=UncondSurvProbDiscrete(clow[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
			}
			
          
			if (zeroL>0)
			
				{  
			    	clow[k]=-(clow[k]);
     			    zeroL=UncondSurvProbDiscrete(clow[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k]; 
					
					while (zeroL>0)
					{
						alpha=1.2*alpha;
						clow[k]=(clow[k])*alpha;;
						zeroL=UncondSurvProbDiscrete(clow[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k]; 
					}
				}
	  
		cmid[k]=(clow[k]+chigh[k])/2;

		zero=UncondSurvProbDiscrete(cmid[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];	  
	  }



   else
		
		{		  		  
		  clow[k]=ci[k]*alpha;
		  chigh[k]=ci[k]/alpha;
		  cmid[k]=(clow[k]+chigh[k])/2;

    	  double zeroL=UncondSurvProbDiscrete(clow[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k]; 		  

        	while (zeroL>0)

				{
    				alpha=1.2*alpha;
					clow[k]=clow[k]*alpha;
					zeroL=UncondSurvProbDiscrete(clow[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k]; 		  			
				}
			
		  double zeroH=UncondSurvProbDiscrete(chigh[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
	
		  alpha=1.5; 

				  
			if (zeroH<0)

			{
    	    	alpha=1.2*alpha;
				chigh[k]=chigh[k]/alpha;
				zeroH=UncondSurvProbDiscrete(chigh[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
			}
			
          
			if (zeroH<0)
			
				{  
			    	chigh[k]=-(chigh[k]);
	
					 zeroH=UncondSurvProbDiscrete(chigh[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
					
					while (zeroH<0)
					{
					   alpha=1.2*alpha;
						chigh[k]=(chigh[k])*alpha;;
     					zeroH=UncondSurvProbDiscrete(chigh[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  

					}
				}
	  
		
		cmid[k]=(clow[k]+chigh[k])/2;

		zero=UncondSurvProbDiscrete(cmid[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];  
   }

		while (fabs(zero)>0.000001)

		{
	   		 if (zero>0)

				{	
					chigh[k]=cmid[k];
					cmid[k]=(chigh[k]+clow[k])/2;
				}

			else 

				{
					clow[k]=cmid[k];
					cmid[k]=(chigh[k]+clow[k])/2;
				}

			   
			double test	=UncondSurvProbDiscrete(cmid[k],UncondMean,UncondStdDev,probMarket,theta);
			double test2=DefaultProb[k];
			zero=UncondSurvProbDiscrete(cmid[k],UncondMean,UncondStdDev,probMarket,theta)-DefaultProb[k];
	 		
		 }
 
 	ci[k]=cmid[k];

	
  }

			
	Matrix<double> CondLoss(MaxNoUnits+1, NoPoints);
	
	
	Array<double> LossTemp(MaxNoUnits+1);

	double dens;
	
	
	for (int i=0; i<NoPoints; i++)

	{
	  
	   Array<double> condDefaultProbi(NoNames);
	   

	   for (int j=0; j<NoNames; j++)

		{
		
    	  condDefaultProbi[j]=GtoNormalCum((ci[j]-theta[i]-UncondMean)/UncondStdDev);				
			
		}
		 
		
		dens=probMarket[i];
			   
		LossTemp = dens*(CondLossDist(NoNames, LossGivenDefault, condDefaultProbi));

		
		for (int k=0; k<=MaxNoUnits; k++)

		{ 
			CondLoss[k][i] = LossTemp[k];
		}	
	
	
	}


	
	for (j=0; j<=MaxNoUnits; j++)
	{
		
	   for (int i=0; i<(NoPoints); i++)

		{ 		
		    Loss[j] = Loss[j]+CondLoss[j][i];
		}
			
		
	}

	
	
	return Loss;
}


// Calc distribution in basic random factor regime where X = beta(i)I(Theta< z < Theta)

Array<double> LossDistRandomNFactorsSimp(
	 const int	   NoNames,
	 Array<double>  beta,
	 Array<double> theta,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb,
	 double		    NoIntervals,
	 double			lbound,
	 double			ubound,
	 double			mbound)
{ 
	
	Array<int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
	int MaxNoUnits=LossParam[1];
	
	Array<double> Loss(MaxNoUnits+1);	  
	Array<double> ci(NoNames);
	
	Array<double> LossProb(MaxNoUnits+1);

	int NoBeta=beta.size();

    // calc mean and variance
   	double UncondMean=beta[0]*GtoNormalDen(theta[0]);
	double V=beta[0]*beta[0]*(GtoNormalCum(theta[0])-theta[0]*GtoNormalDen(theta[0]));	 
	for (int s=1; s<(NoBeta-1); s++)
	{
		UncondMean=UncondMean+beta[s]*(GtoNormalDen(theta[s])-GtoNormalDen(theta[s-1]));
		V=V+beta[s]*beta[s]*(GtoNormalCum(theta[s])-GtoNormalCum(theta[s-1])+theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*GtoNormalDen(theta[s]));
	}
	UncondMean=UncondMean+beta[NoBeta-1]*(-GtoNormalDen(theta[NoBeta-2]));
	V=V+beta[NoBeta-1]*beta[NoBeta-1]*(1-GtoNormalCum(theta[NoBeta-2])+theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2])); 
    V=V-UncondMean*UncondMean;
    // fail if variance > 1
	if (V>1)
    {
        for (int j=0; j<MaxNoUnits+1; j++)
            Loss[j]=-1;
        return Loss;
    }
  	double UncondStdDev=sqrt(1-V);


    // calibrate individual name loss threshold
	for (int k=0; k<NoNames; k++)
        GtoNormalCumInv(DefaultProb[k],&ci[k]);

	Array<double> clow(NoNames);
	Array<double> chigh(NoNames);
	Array<double> cmid(NoNames);
    double alpha=1.05;
    double zero;

    for (k=0; k<NoNames; k++)
	{ 
        // bracket c first
        if (ci[k]>0)
        {
            clow[k]=ci[k]/alpha;
            chigh[k]=ci[k]*alpha;
            cmid[k]=(clow[k]+chigh[k])/2;

            // locate upper bound
            double zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
        	while (zeroH<0)
            {
                alpha=1.2*alpha;
                chigh[k]=chigh[k]*alpha;
                zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
            }

            // locate lower bound			
            alpha = 1.05;
        	double zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	
			if (zeroL>0) // what is this? some sort of heuristics?
			{
    	    	alpha=1.2*alpha;
				clow[k]=clow[k]/alpha;
				zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
			}
			if (zeroL>0)
			{  
                clow[k]=-(clow[k]);
                zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 	
				while (zeroL>0)
				{
					alpha=1.2*alpha;
					clow[k]=(clow[k])*alpha;;
					zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 
				}
			}
	  
    		cmid[k]=(clow[k]+chigh[k])/2;
            zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	  
        } else {
            clow[k]=ci[k]*alpha;
            chigh[k]=ci[k]/alpha;
            cmid[k]=(clow[k]+chigh[k])/2;

            double zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  
        	while (zeroL>0)
			{
   				alpha=1.2*alpha;
				clow[k]=clow[k]*alpha;
				zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  			
			}
			
            double zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  	
            alpha=1.5;  
			if (zeroH<0)
			{
    	    	alpha=1.2*alpha;
				chigh[k]=chigh[k]/alpha;
				zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
			}
			if (zeroH<0)
			{  
                chigh[k]=-(chigh[k]);
                zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	
				while (zeroH<0)
				{
                    alpha=1.2*alpha;
                    chigh[k]=(chigh[k])*alpha;;
    			    zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
				}
			}	  
		
		    cmid[k]=(clow[k]+chigh[k])/2;
            zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
        }

        // do bisection
		while (fabs(zero)>0.000001)
		{
            if (zero>0)
            {	
                chigh[k]=cmid[k];
                cmid[k]=(chigh[k]+clow[k])/2;
            } else {
                clow[k]=cmid[k];
				cmid[k]=(chigh[k]+clow[k])/2;
			}
            zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
        }
 
        // found c
 	    ci[k]=cmid[k];	
    }


    // prepare integration boundary
	double lengthL=2*(mbound-lbound)/NoIntervals;
	double lengthH=2*(ubound-mbound)/NoIntervals;
	int ibound=NoIntervals/2;
			
	Matrix<double> CondLoss(MaxNoUnits+1, NoIntervals+1);
	Array<double> LossTemp(MaxNoUnits+1);

	double dens;
	double xi; // market factor
	double betai=beta[0]; // beta under current market factor

    // loop over all integration points and calc integrand value for all points
	for (int i=0; i<NoIntervals+1; i++)
	{
        // calculate current market factor
		if (i<=ibound)
			xi=lbound+i*lengthL;
		else
			xi=lbound+(NoIntervals/2)*lengthL +lengthH*(i-ibound);
	  
        // calc conditional prob for each name under current market factor
        Array<double> condDefaultProbi(NoNames);
        for (int j=0; j<NoNames; j++)
		{
            // locate current beta
			for (int s=0; s<NoBeta; s++)
			{ 
				if ((s==0) && (xi <= theta[0]))
                    betai=beta[0];
				
				if ((s==(NoBeta-1)) && (xi>theta[NoBeta-2]))
                    betai=beta[NoBeta-1];
				
				if ((s>0) && (xi>theta[s-1]) && (xi<=theta[s]))
                    betai=beta[s];
			}
			condDefaultProbi[j]=GtoNormalCum((ci[j]-betai*xi-UncondMean)/UncondStdDev);
		}
		 
		dens=exp(-1*(xi)*(xi)/2)/sqrt(2*3.141592653);			   
		LossTemp = dens*(CondLossDist(NoNames, LossGivenDefault, condDefaultProbi));
		
		for (int k=0; k<=MaxNoUnits; k++)
		{ 
			CondLoss[k][i] = LossTemp[k];
		}	
	}

    // do simpson integration
	double wi=1;
	for (int j=0; j<=MaxNoUnits; j++)
	{
		Loss[j] = 0;
		for (int i=0; i<(NoIntervals+1); i++)
		{ 		
			wi=4/3;

			if ( (i==0) || (i==NoIntervals) ) 
			{wi=27/12;}

			if ( (i==1) || (i==(NoIntervals-1)) ) 
			{wi=0;}

			if ( (i==2) || (i==(NoIntervals-2)) ) 
			{wi=13/12;}

            // is this right, use 4/3 for all points in the middle?
			if (i<=ibound)
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthL;
			else
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthH;
		}
	}
	
	return Loss;
}

Matrix<double> LossDistRandomNFactorsSimpsBuck(

	 const int	   NoNames,
	 double		   fracG1,
	 double		   fracG2,
	 double		   probG1,
	 double		   probG2,
	 double		   probJ,
	 Array<double> beta,
	 Array<double> theta,
	 Array<int>    &Notional,
	 Array<double> &Recovery,
	 Array<double> &DefaultProb,
	 double		    NoIntervals,
	 double			lbound,
	 double			ubound,
	 double			mbound,
	 double			NoUnits)


{ 
	Array<double> LossGivenDef(NoNames);

	int cumLoss=0;
	
	int NoLossUnits=ceil(NoUnits)+1;

	for (int j=0; j<NoNames; j++)
	{
	  LossGivenDef[j]=Notional[j]*(1-Recovery[j]);
	  cumLoss=cumLoss+LossGivenDef[j];
	}

	
	double LossUnit=cumLoss/NoUnits;

	for (j=0; j<NoNames; j++)
	{
	  LossGivenDef[j]=LossGivenDef[j]/LossUnit;
	}
	
	

	Matrix<double> LossPr(NoLossUnits+1,2);	 
	
	Array<double> Mean(NoLossUnits+1);	 

	Array<double> ci(NoNames);
	
	Array<double> LossProb(NoLossUnits+1);

	int NoBeta=beta.size();
	
	double z;

	double lengthL=2*(mbound-lbound)/NoIntervals;

	double lengthH=2*(ubound-mbound)/NoIntervals;

	int ibound=NoIntervals/2;

	double NoNamesG1d=fracG1*NoNames;
	double NoNamesG2d=fracG2*NoNames;
	
	int NoNamesG1 = (int)NoNamesG1d;
	int NoNamesG2 = (int)NoNamesG2d;


	for (int k=0; k<NoNames; k++)

	{ 
		if (k<NoNamesG1)     
		{
			DefaultProb[k]=(DefaultProb[k]-probG1-probJ)/(1-probG1-probJ);
		}
		
		else
		{
		 if (k<(NoNamesG2+NoNamesG1))     
			{
			DefaultProb[k]=(DefaultProb[k]-probG2-probJ)/(1-probG2-probJ);
			}
		 else
			{
			DefaultProb[k]=(DefaultProb[k]-probJ)/(1-probJ);
			}

		}
	}
  


	for (k=0; k<NoNames; k++)

	{ 
      GtoNormalCumInv(DefaultProb[k],&z);
	  ci[k]=z;
	}
  
	
	double UncondMean=beta[0]*GtoNormalDen(theta[0]);
	double V=beta[0]*beta[0]*(GtoNormalCum(theta[0])-theta[0]*GtoNormalDen(theta[0]));
	 
	for (int s=1; s<(NoBeta-1); s++)

	{
		UncondMean=UncondMean+beta[s]*(GtoNormalDen(theta[s])-GtoNormalDen(theta[s-1]));
		V=V+beta[s]*beta[s]*(GtoNormalCum(theta[s])-GtoNormalCum(theta[s-1])+theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*GtoNormalDen(theta[s]));
	}
	
		UncondMean=UncondMean+beta[NoBeta-1]*(-GtoNormalDen(theta[NoBeta-2]));
		V=V+beta[NoBeta-1]*beta[NoBeta-1]*(1-GtoNormalCum(theta[NoBeta-2])+theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2])); 
	
	
		 V=V-UncondMean*UncondMean;

		 if (V>1)
		 
		 {
			 for (int j=0; j<(NoLossUnits+1); j++)
			 {
				 LossPr[j][0]=-1;
				 LossPr[j][1]=-1;
			 }
		 return LossPr;
		 }


	double UncondStdDev=sqrt(1-V);

	Array<double> clow(NoNames);
	Array<double> chigh(NoNames);
	Array<double> cmid(NoNames);
	double alpha=1.05;
	double zero;
	

for (k=0; k<NoNames; k++)

	{ 
	  if (ci[k]>0)
		
	  {
		clow[k]=ci[k]/alpha;
		chigh[k]=ci[k]*alpha;
		cmid[k]=(clow[k]+chigh[k])/2;

		double zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  

        	while (zeroH<0)

				{
    				alpha=1.2*alpha;
					chigh[k]=chigh[k]*alpha;
					zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
		  		}
			
    	double zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
	
    	  alpha=1.05; 

				  
			if (zeroL>0)

			{
    	    	alpha=1.2*alpha;
				clow[k]=clow[k]/alpha;
				zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
			}
			
          
			if (zeroL>0)
			
				{  
			    	clow[k]=-(clow[k]);
     			    zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 
					
					while (zeroL>0)
					{
						alpha=1.2*alpha;
						clow[k]=(clow[k])*alpha;;
						zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 
					}
				}
	  
		cmid[k]=(clow[k]+chigh[k])/2;

		zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	  
	  }



   else
		
		{		  		  
		  clow[k]=ci[k]*alpha;
		  chigh[k]=ci[k]/alpha;
		  cmid[k]=(clow[k]+chigh[k])/2;

    	  double zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  

        	while (zeroL>0)

				{
    				alpha=1.2*alpha;
					clow[k]=clow[k]*alpha;
					zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  			
				}
			
		  double zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
	
		  alpha=1.5; 

				  
			if (zeroH<0)

			{
    	    	alpha=1.2*alpha;
				chigh[k]=chigh[k]/alpha;
				zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
			}
			
          
			if (zeroH<0)
			
				{  
			    	chigh[k]=-(chigh[k]);
	
					 zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
					
					while (zeroH<0)
					{
					   alpha=1.2*alpha;
						chigh[k]=(chigh[k])*alpha;;
     					zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  

					}
				}
	  
		
		cmid[k]=(clow[k]+chigh[k])/2;

		zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
   }

		while (fabs(zero)>0.000001)

		{
	   		 if (zero>0)

				{	
					chigh[k]=cmid[k];
					cmid[k]=(chigh[k]+clow[k])/2;
				}

			else 

				{
					clow[k]=cmid[k];
					cmid[k]=(chigh[k]+clow[k])/2;
				}

			   zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
	 		
		 }
 
 	ci[k]=cmid[k];


  }

			
	Matrix<double> CondLoss(NoLossUnits+1, NoIntervals+1);
	
	Matrix<double> CondExpected(NoLossUnits+1, NoIntervals+1);
	

	Matrix<double> LossTemp(NoLossUnits+1,2);

	double dens;
	double xi;
	double betai=beta[0];

	for (int i=0; i<NoIntervals+1; i++)

	{
		if (i<=ibound)
		{
			xi=lbound+i*lengthL;
		}

		else
		{ 
			xi=lbound+(NoIntervals/2)*lengthL +lengthH*(i-ibound);
		}

	  
	   Array<double> condDefaultProbi(NoNames);
	   

	   for (int j=0; j<NoNames; j++)

		{
		
			for (int s=0; s<NoBeta; s++)

			{ 
				if ((s==0) && (xi <= theta[0]))
				{ 
				  betai=beta[0];
				}
				
				if ((s==(NoBeta-1)) && (xi>theta[NoBeta-2]))
				
				{		 
				  betai=beta[NoBeta-1];
				}
				
				if ((s>0) && (xi>theta[s-1]) && (xi<=theta[s]))
				{
				  betai=beta[s];
				}
			}


			condDefaultProbi[j]=GtoNormalCum((ci[j]-betai*xi-UncondMean)/UncondStdDev);				
			
		}
		 
		
		dens=exp(-1*(xi)*(xi)/2)/sqrt(2*3.141592653);
			   
		LossTemp = (CondLossDistUnits(NoLossUnits, LossGivenDef, condDefaultProbi));


		
		for (int k=0; k<=NoLossUnits; k++)

		{ 
			CondLoss[k][i] = (1-probG1-probG2-probJ)*LossTemp[k][0]*dens;
			CondExpected[k][i] =  (1-probG1-probG2-probJ)*LossTemp[k][0]*LossTemp[k][1]*dens;
		}	
	
	
	}



	double wi=1;

	for (j=0; j<=NoLossUnits; j++)
	{
		
		LossPr[j][0] = 0;

		for (int i=0; i<(NoIntervals+1); i++)

		{ 		
			wi=4/3;

			if ( (i==0) || (i==NoIntervals) ) 
			{wi=27/12;}

			if ( (i==1) || (i==(NoIntervals-1)) ) 
			{wi=0;}

			if ( (i==2) || (i==(NoIntervals-2)) ) 
			{wi=13/12;}


			if (i<=ibound)	
			{
				LossPr[j][0] = LossPr[j][0]+wi*CondLoss[j][i]*lengthL;
				LossPr[j][1]=LossPr[j][1]+wi*CondExpected[j][i]*lengthL;
			}

			else
			{
				LossPr[j][0] = LossPr[j][0]+wi*CondLoss[j][i]*lengthH;
				LossPr[j][1]  =LossPr[j][1]+wi*CondExpected[j][i]*lengthH;
			}
			
		}

		if (LossPr[j][0]==0)
		{ 
			LossPr[j][1]=0.5*(2*j-1);
		}

		else
		{
		LossPr[j][1]=LossPr[j][1]/LossPr[j][0];
		}
	
	}

	//Loss[NoNames]=Loss[NoNames]+probJ;
	

	return LossPr;
}
















Matrix<double> LossDistRandomFactorGeneralMultper(

	 const int	   NoNames,
	 double		   beta,
	 double		   DispGen,
	 double		   Disp1,
	 double		   Disp2,
	 double		   theta1,
	 double		   theta2,
	 Array<int>    &LossGivenDefault,
	 Matrix<double> &DefaultProb)


{ 
	Array<int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
	
	int MaxNoUnits=LossParam[1];
	
	int NoPeriods = DefaultProb.colSize();

	Matrix<double> Loss(MaxNoUnits+1,NoPeriods);	  
	Matrix<double> ci;

	ci.resize(NoNames,NoPeriods);

	Array<double> LossProb(MaxNoUnits+1);
	
	double z;

	for (int j=0; j<NoPeriods; j++)
	{
		for (int v=0; v<NoNames; v++)

		{ 
			GtoNormalCumInv(DefaultProb[v][j],&z);
			ci[v][j]=z;
		}
	}

	double weights[40];
	double points[40];

	imsl_d_gauss_quad_rule(40,weights, points, IMSL_HERMITE,0);

	
	if (theta1>theta2)

	{
	  double temp=theta2;
	  theta2=theta1;
	  theta1=temp;

	  double temp2=Disp2;
	  Disp2=Disp1;
	  Disp1=temp2;

	}
		
	 
	 double UncondMean = 2*Disp2*GtoNormalDen(theta2)+2*Disp1*GtoNormalDen(theta1);
	 double V;
	 
	 V=beta*beta+DispGen*DispGen+(Disp1+Disp2)*(Disp1+Disp2)-2*beta*(Disp1+Disp2);
	 V=V+4*Disp1*(Disp2+beta)*(-theta1*GtoNormalDen(theta1)+GtoNormalCum(theta1))+4*Disp2*(Disp1-beta)*(theta2*GtoNormalDen(theta2)-GtoNormalCum(theta2));
	 V=V-UncondMean*UncondMean;

	 if (V>1)
		 
	 {
		 for (int u=0; u<NoPeriods; u++)
		 {
			for (int j=0; j<MaxNoUnits+1; j++)
			 {
				Loss[j][u]=-1;
			}
		 }
	 
	 return Loss;
	 
	 }


	 double UncondStdDev=sqrt(1-V);


	Matrix<double> clow(NoNames,NoPeriods);
	Matrix<double> chigh(NoNames,NoPeriods);
	Matrix<double> cmid(NoNames,NoPeriods);
	
	double alpha=1.2;
	

	for (j=0; j<NoPeriods; j++)

	{

		for (int k=0; k<NoNames; k++)

		{ 
			if (ci[k][j]>0)
		
			{
				clow[k][j]=ci[k][j]/alpha;
				chigh[k][j]=ci[k][j]*alpha;
				cmid[k][j]=(clow[k][j]+chigh[k][j])/2;

				double zero=UncondSurvProbGen(cmid[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];

				while (fabs(zero)>0.000001)
				{
	   	 
				if (zero>0)

					{	
						chigh[k][j]=cmid[k][j];
						cmid[k][j]=(chigh[k][j]+clow[k][j])/2;
					}

				else 

					{
						clow[k][j]=cmid[k][j];
						cmid[k][j]=(chigh[k][j]+clow[k][j])/2;
					}
	
				 zero=UncondSurvProbGen(cmid[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];
	
				}
			}

			else
		
			{
				clow[k][j]=ci[k][j]*alpha;
				chigh[k][j]=ci[k][j]/alpha;
	   
				double zeroL=UncondSurvProbGen(clow[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];  
		  
				while (zeroL>0)

				{
    				alpha=1.2*alpha;
					clow[k][j]=clow[k][j]*alpha;
					zeroL=UncondSurvProbGen(clow[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];  
		  
				}


				double zeroH=UncondSurvProbGen(chigh[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];  
	
				alpha=1.5; 
		  
				if (zeroH<0)

				{
    				alpha=1.2*alpha;
					chigh[k][j]=chigh[k][j]/alpha;
					zeroH=UncondSurvProbGen(chigh[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];    
				}
			
          
				if (zeroH<0)
			
				{  
					chigh[k][j]=-(chigh[k][j]);
		
					while (zeroH<0)
					{
						alpha=1.2*alpha;
						chigh[k][j]=(chigh[k][j])*alpha;;
						zeroH=UncondSurvProbGen(chigh[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];  
					 }
				}
	  
			}

			cmid[k][j]=(clow[k][j]+chigh[k][j])/2;


			double zero=UncondSurvProbGen(cmid[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];
	
	    
			while (fabs(zero)>0.000001)
	
			{
	   		 if (zero>0)

				{	
					chigh[k][j]=cmid[k][j];
					cmid[k][j]=(chigh[k][j]+clow[k][j])/2;
				}

			 else 

				{
					clow[k][j]=cmid[k][j];
					cmid[k][j]=(chigh[k][j]+clow[k][j])/2;
					double a=cmid[k][j];
					double a1=chigh[k][j];
					double a2=clow[k][j];
				}

			  zero=UncondSurvProbGen(cmid[k][j],UncondMean,UncondStdDev,beta,DispGen,Disp1,Disp2,theta1,theta2)-DefaultProb[k][j];
			}
	  
		 
 		ci[k][j]=cmid[k][j];
	
		}
	}
	

	for (int l=0; l<NoPeriods; l++)

	{

	Matrix<double> CondLoss(MaxNoUnits+1, 40);
	Array<double> LossTemp(MaxNoUnits+1);

	for (int i=0; i<40; i++)

	{	
		Array<double> condDefaultProbi(NoNames);
		Array<double> condDefaultProbi2(NoNames);

		for (int j=0; j<NoNames; j++)

		{
			if ((points[i]*sqrt(2))>theta2)
			{
			condDefaultProbi[j]=GtoNormalCum((ci[j][l]-(beta+DispGen-Disp1-Disp2)*points[i]*sqrt(2)-UncondMean)/UncondStdDev);	
     		condDefaultProbi2[j]=GtoNormalCum((ci[j][l]-(beta-DispGen-Disp1-Disp2)*points[i]*sqrt(2)-UncondMean)/UncondStdDev);
			}
			else
			{
				if ((points[i]*sqrt(2))>theta1)
				{ 
				condDefaultProbi[j]=GtoNormalCum((ci[j][l]-(beta+DispGen-Disp1+Disp2)*points[i]*sqrt(2)-UncondMean)/UncondStdDev);	
				condDefaultProbi2[j]=GtoNormalCum((ci[j][l]-(beta-DispGen-Disp1+Disp2)*points[i]*sqrt(2)-UncondMean)/UncondStdDev);
				}
				else 
				{
				condDefaultProbi[j]=GtoNormalCum((ci[j][l]-(beta+DispGen+Disp1+Disp2)*points[i]*sqrt(2)-UncondMean)/UncondStdDev);
				condDefaultProbi2[j]=GtoNormalCum((ci[j][l]-(beta-DispGen+Disp1+Disp2)*points[i]*sqrt(2)-UncondMean)/UncondStdDev);
				}		
			}
		}
		 
		LossTemp = 0.5*(CondLossDist(NoNames, LossGivenDefault, condDefaultProbi)+CondLossDist(NoNames, LossGivenDefault, condDefaultProbi2));

		for (int k=0; k<=MaxNoUnits; k++)

		{ 
			CondLoss[k][i] = LossTemp[k];
		}	
		
	}

	for (int j=0; j<=MaxNoUnits; j++)
	{
		
		Loss[j][l] = 0;

		for (int i=0; i<40; i++)

		{ 
			Loss[j][l] = Loss[j][l]+weights[i]*CondLoss[j][i]/sqrt(3.141592653);
		}
	}
	
	}
	
	return Loss;

}





Array<double> UncondMomentsNfactor(

	 Array<double>  beta,
	 Array<double>  theta)


{

	int NoBeta=beta.size();

	Array<double> moments(3);  
	
	 double UncondMean=beta[0]*GtoNormalDen(theta[0]);
	 double V=beta[0]*beta[0]*(GtoNormalCum(theta[0])-theta[0]*GtoNormalDen(theta[0]));
	 
	 double Sk=beta[0]*beta[0]*beta[0]*(-2*GtoNormalDen(theta[0])-theta[0]*theta[0]*GtoNormalDen(theta[0]));	

	 for (int s=1; s<(NoBeta-1); s++)

	{
		UncondMean=UncondMean+beta[s]*(GtoNormalDen(theta[s])-GtoNormalDen(theta[s-1]));
		V=V+beta[s]*beta[s]*(GtoNormalCum(theta[s])-GtoNormalCum(theta[s-1])+theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*GtoNormalDen(theta[s]));
		Sk=Sk+beta[s]*beta[s]*beta[s]*(2*GtoNormalDen(theta[s-1])-2*GtoNormalDen(theta[s])+theta[s-1]*theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*theta[s]*GtoNormalDen(theta[s]));	
	}
	
		UncondMean=UncondMean+beta[NoBeta-1]*(-GtoNormalDen(theta[NoBeta-2]));
		V=V+beta[NoBeta-1]*beta[NoBeta-1]*(1-GtoNormalCum(theta[NoBeta-2])+theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2])); 
		Sk=Sk+beta[NoBeta-1]*beta[NoBeta-1]*beta[NoBeta-1]*(2*GtoNormalDen(theta[NoBeta-2])+theta[NoBeta-2]*theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2]));	
	
	moments[0]=UncondMean;
	moments[1]=V-UncondMean*UncondMean;
	moments[2]=Sk-3*UncondMean*V-2*UncondMean*UncondMean*UncondMean;
	
	return moments;


}





Array<double> LossDistRandomNFactorsHomogenSimp(

	 const int	    NoNames,
	 Array<double>  beta,
	 Array<double>  theta,
	 double			DefaultProbT,
	 double		    NoIntervals,
	 double			lbound,
	 double			ubound,
	 double			mbound)


{

	int MaxNoUnits=NoNames;
	
	Array<double> Loss(MaxNoUnits+1);	  
	double ci;
	
	Array<double> LossProb(MaxNoUnits+1);

	int NoBeta=beta.size();
	
	double z;

	double lengthL=2*(mbound-lbound)/NoIntervals;

	double lengthH=2*(ubound-mbound)/NoIntervals;

	int ibound=NoIntervals/2;

	
	double DefaultProb=DefaultProbT;

     GtoNormalCumInv(DefaultProb,&z);
	 ci=z;
	  
	
	 double UncondMean=beta[0]*GtoNormalDen(theta[0]);
	  
	 double V=beta[0]*beta[0]*(GtoNormalCum(theta[0])-theta[0]*GtoNormalDen(theta[0]));
	 
	for (int s=1; s<(NoBeta-1); s++)

	{
		UncondMean=UncondMean+beta[s]*(GtoNormalDen(theta[s])-GtoNormalDen(theta[s-1]));
		V=V+beta[s]*beta[s]*(GtoNormalCum(theta[s])-GtoNormalCum(theta[s-1])+theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*GtoNormalDen(theta[s]));
	}
	
		UncondMean=UncondMean+beta[NoBeta-1]*(-GtoNormalDen(theta[NoBeta-2]));
		V=V+beta[NoBeta-1]*beta[NoBeta-1]*(1-GtoNormalCum(theta[NoBeta-2])+theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2])); 
	
	
		V=V-UncondMean*UncondMean;

		 
		 if (V>1)
		 
		 {
			 for (int j=0; j<MaxNoUnits+1; j++)
			 {
				 Loss[j]=-1;
			 }
		 return Loss;
		 }

	double UncondStdDev=sqrt(1-V);

	double clow;
	double chigh;
	double cmid;
    double alpha=1.1;
    double zero;
	

  if (ci>0)
		
	  {
		clow=ci/alpha;
		chigh=ci*alpha;
		cmid=(clow+chigh)/2;

		double zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,beta,theta)-DefaultProb; 		  

        	while (zeroH<0)

				{
    				alpha=1.2*alpha;
					chigh=chigh*alpha;
					zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
		  		}
			
    	double zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
	
    	  alpha=1.1; 

				  
			if (zeroL>0)

			{
    	    	alpha=1.2*alpha;
				clow=clow/alpha;
				zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
			}
			
          
			if (zeroL>0)
			
				{  
			    	clow=-(clow);
					zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
					
					while (zeroL>0)
					{
						alpha=1.2*alpha;
						clow=(clow)*alpha;;
						zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
					}
				}
	  
		cmid=(clow+chigh)/2;

		zero=UncondSurvProbNfact(cmid,UncondMean,UncondStdDev,beta,theta)-DefaultProb;	  
	  }



   else
		
		{		  		  
		  clow=ci*alpha;
		  chigh=ci/alpha;
		  cmid=(clow+chigh)/2;

    	  double zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,beta,theta)-DefaultProb; 		  

        	while (zeroL>0)

				{
    				alpha=1.2*alpha;
					clow=clow*alpha;
					zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,beta,theta)-DefaultProb; 		  			
				}
			
		  double zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
	
		  alpha=1.5; 

				  
			if (zeroH<0)

			{
    	    	alpha=1.2*alpha;
				chigh=chigh/alpha;
				zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
			}
			
          
			if (zeroH<0)
			
				{  
			    	chigh=-(chigh);
	
					 zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
					
					while (zeroH<0)
					{
					   alpha=1.2*alpha;
						chigh=(chigh)*alpha;;
     					zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  

					}
				}
	  
		
		cmid=(clow+chigh)/2;

		zero=UncondSurvProbNfact(cmid,UncondMean,UncondStdDev,beta,theta)-DefaultProb;  
   }

		while (fabs(zero)>0.000001)

		{
	   		 if (zero>0)

				{	
					chigh=cmid;
					cmid=(chigh+clow)/2;
				}

			else 

				{
					clow=cmid;
					cmid=(chigh+clow)/2;
				}

			   zero=UncondSurvProbNfact(cmid,UncondMean,UncondStdDev,beta,theta)-DefaultProb;
	 		
		 }
 
 	ci=cmid;

	
			
	Matrix<double> CondLoss(MaxNoUnits+1, NoIntervals+1);
	
	
	Array<double> LossTemp(MaxNoUnits+1);

	
	double dens;
	double xi;
	double betai=beta[0];

	for (int i=0; i<NoIntervals+1; i++)

	{
		if (i<=ibound)
		{
			xi=lbound+i*lengthL;
		}

		else
		{ 
			xi=lbound+(NoIntervals/2)*lengthL +lengthH*(i-ibound);
		}

	  
	   
	   double condDefaultProbi;
	   


    	for (int s=0; s<NoBeta; s++)

			{ 
				if ((s==0) && (xi <= theta[0]))
				{ 
				  betai=beta[0];
				}
				
				if ((s==(NoBeta-1)) && (xi>theta[NoBeta-2]))
				
				{		 
				  betai=beta[NoBeta-1];
				}
				
				if ((s>0) && (xi>theta[s-1]) && (xi<=theta[s]))
				{
				  betai=beta[s];
				}
			}


		condDefaultProbi=GtoNormalCum((ci-betai*xi-UncondMean)/UncondStdDev);				
			
		dens=exp(-1*(xi)*(xi)/2)/sqrt(2*3.141592653);
			   
		
		for (int k=0; k<MaxNoUnits; k++)

		
		{
		  	
			double comb=imsl_d_gamma(MaxNoUnits+1)/imsl_d_gamma(k+1)/imsl_d_gamma(MaxNoUnits-k+1);
			
			LossTemp[k]= dens*comb*pow(condDefaultProbi,k)*pow(1-condDefaultProbi,MaxNoUnits-k);

		}		
		
		
		double comb=imsl_d_gamma(MaxNoUnits+1)/imsl_d_gamma(MaxNoUnits+1)/imsl_d_gamma(1);
			
		LossTemp[MaxNoUnits]= dens*(comb*pow(condDefaultProbi,MaxNoUnits)*pow(1-condDefaultProbi,0));

		double test=imsl_d_gamma(3);
		
		for ( k=0; k<=MaxNoUnits; k++)

		{ 
			CondLoss[k][i] = LossTemp[k];
		}	
	
	
	}

	double wi=1;

	for (int j=0; j<=MaxNoUnits; j++)
	{
		
		Loss[j] = 0;

		for (int i=0; i<(NoIntervals+1); i++)

		{ 		
			wi=4/3;

			if ( (i==0) || (i==NoIntervals) ) 
			{wi=27/12;}

			if ( (i==1) || (i==(NoIntervals-1)) ) 
			{wi=0;}

			if ( (i==2) || (i==(NoIntervals-2)) ) 
			{wi=13/12;}


			if (i<=ibound)	
			{
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthL;
			}

			else
			{
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthH;
			}
			
		}
	}


	return Loss;

}


// Calc loss distribution for homogenous portfolio uing RR

Array<double> LossDistRandomNFactorsHomogenSimpRR(
	 const int	    NoNames,
	 Array<double>  Beta,
	 Array<double>  Theta,
	 double			ProbJ,
	 double			DefaultProbT,
	 double		    NoIntervals,
	 Array<double>	Recovery,
	 Array<double>	&ThresholdR,
	 double			lbound,
	 double			ubound,
	 double			mbound)
{
	double loss0=100*(1-Recovery[0]);
	double loss1=100*(1-Recovery[1]);
	double loss2=100*(1-Recovery[2]);

	int minloss=GCD(int(loss0),GCD(int(loss1),int(loss2)));
	int rat0=loss0/minloss;
	int rat1=loss1/minloss;
	int rat2=loss2/minloss;
	int MaxNoUnits=NoNames*rat0;

    Array<double> Loss(MaxNoUnits+1);

    // sanity check
    int i;
    for (i = 0; i < Theta.size() - 1; ++i) {
        if (Theta[i+1] < Theta[i]) {
            for (int j=0; j<MaxNoUnits+1; j++)
            {
                Loss[j]=-1;
            }
            return Loss;
        }
    }
    for (i = 0; i < ThresholdR.size() - 1; ++i) {
        if (ThresholdR[i+1] < ThresholdR[i]) {
            for (int j=0; j<MaxNoUnits+1; j++)
            {
                Loss[j]=-1;
            }
            return Loss;
        }
    }

	int NoBeta=Beta.size();
	int NoRec=Recovery.size();
	
    // Calc mean and variance

    double UncondMean=Beta[0]*GtoNormalDen(Theta[0]);
    double V=Beta[0]*Beta[0]*(GtoNormalCum(Theta[0])-Theta[0]*GtoNormalDen(Theta[0]));
    for (int s=1; s<(NoBeta-1); s++)
	{
		UncondMean=UncondMean+Beta[s]*(GtoNormalDen(Theta[s])-GtoNormalDen(Theta[s-1]));
		V=V+Beta[s]*Beta[s]*(GtoNormalCum(Theta[s])-GtoNormalCum(Theta[s-1])+Theta[s-1]*GtoNormalDen(Theta[s-1])-Theta[s]*GtoNormalDen(Theta[s]));
	}
	UncondMean=UncondMean+Beta[NoBeta-1]*(-GtoNormalDen(Theta[NoBeta-2]));
	V=V+Beta[NoBeta-1]*Beta[NoBeta-1]*(1-GtoNormalCum(Theta[NoBeta-2])+Theta[NoBeta-2]*GtoNormalDen(Theta[NoBeta-2])); 
    V=V-UncondMean*UncondMean;	 
    if (V>1)
    {
        for (int j=0; j<MaxNoUnits+1; j++)
        {
            Loss[j]=-1;
        }
        return Loss;
    }
	double UncondStdDev=sqrt(1-V);

    // Calc default threshold ci, what is this ProJ thing?

    double DefaultProb=(DefaultProbT-ProbJ)/(1-ProbJ);
    double ci;
    GtoNormalCumInv(DefaultProb,&ci);

	double clow;
	double chigh;
	double cmid;
    double alpha=1.1;
    double zero;
	
    // bracket c first
    if (ci>0)
    {
		clow=ci/alpha;
		chigh=ci*alpha;
		cmid=(clow+chigh)/2;

		double zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;
       	while (zeroH<0)
		{
			alpha=1.2*alpha;
			chigh=chigh*alpha;
			zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
  		}			
	   	double zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
        alpha=1.1;
        if (zeroL>0)
		{
   	    	alpha=1.2*alpha;
			clow=clow/alpha;
			zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
		}
        if (zeroL>0)
		{  
		   	clow=-(clow);
			zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
			while (zeroL>0)
			{
				alpha=1.2*alpha;
				clow=(clow)*alpha;;
				zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
			}
		}
		cmid=(clow+chigh)/2;
		zero=UncondSurvProbNfact(cmid,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;	  
    } else {
        clow=ci*alpha;
        chigh=ci/alpha;
        cmid=(clow+chigh)/2;

        double zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;
        while (zeroL>0)
		{
			alpha=1.2*alpha;
			clow=clow*alpha;
			zeroL=UncondSurvProbNfact(clow,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb; 		  			
		}
        double zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;
        alpha=1.5;
		if (zeroH<0)
		{
   	    	alpha=1.2*alpha;
			chigh=chigh/alpha;
			zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
		}
        if (zeroH<0)
		{  
	    	chigh=-(chigh);
            zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;
            while (zeroH<0)
			{
                alpha=1.2*alpha;
				chigh=(chigh)*alpha;;
     			zeroH=UncondSurvProbNfact(chigh,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;
			}
		}		
		cmid=(clow+chigh)/2;
		zero=UncondSurvProbNfact(cmid,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;  
    }

    // do bisection
	while (fabs(zero)>0.00001 && (cmid - clow > 0.0001))
	{
        if (zero>0)
		{	
			chigh=cmid;
			cmid=(chigh+clow)/2;
        } else {
			clow=cmid;
			cmid=(chigh+clow)/2;
		}
        zero=UncondSurvProbNfact(cmid,UncondMean,UncondStdDev,Beta,Theta)-DefaultProb;
    }
 
    // found Ci
 	ci=cmid;	
			
	// Calibrate Recovery threshold r2 to match expected recovery Recovery[1]
    
	double tol=0.0001;	
	Matrix<double> Param(NoBeta+NoRec-1,3);
	double thresL=ThresholdR[0]; // assume two thresholds
	double thresH=ThresholdR[1];
	double thresM=(thresL+thresH)/2;
	ThresholdR[1]=thresM; // initial guess
	
	Param=OrderBetaRecTh(Beta,Theta,Recovery,ThresholdR);
	double zeroR=UncondSurvProbNBetaNRec(ci,UncondMean,UncondStdDev,Param)-Recovery[1];
	while ((fabs(zeroR)) >tol && (ThresholdR[1] - ThresholdR[0] > 0.00001))
	{
		if (zeroR<0) // expected recovery too small, reduce r2
		{
			thresH=thresM;
			thresM=(thresL+thresH)/2;	
        } else { // expected recovery too big, increase r2
			thresL=thresM;
			thresM=(thresL+thresH)/2;	
		}	

        ThresholdR[1]=thresM;
        Param=OrderBetaRecTh(Beta,Theta,Recovery,ThresholdR);
        zeroR=UncondSurvProbNBetaNRec(ci,UncondMean,UncondStdDev,Param)-Recovery[1];
	}
    if (fabs(zeroR) > tol) { // can not find an appropriate r2
        for (int j=0; j<MaxNoUnits+1; j++)
        {
            Loss[j]=-1;
        }
        return Loss;
    }

    // dump(ThresholdR, "threshold r");

    // do numerical integration

	Matrix<double> CondLoss(MaxNoUnits+1, NoIntervals+1);
	Array<double> LossTemp(MaxNoUnits+1);

   	double lengthL=2*(mbound-lbound)/NoIntervals;
	double lengthH=2*(ubound-mbound)/NoIntervals;
	int ibound=NoIntervals/2;

	
	double dens;
	double xi;
	double betai=Beta[0];
	double recoveryi=Recovery[0];

	for (i=0; i<NoIntervals+1; i++)
	{
        //i = 350; // xxx

		if (i<=ibound)
			xi=lbound+i*lengthL;
		else
			xi=lbound+(NoIntervals/2)*lengthL +lengthH*(i-ibound);	  
	   
        double condDefaultProbi;
	   
        // locate beta and recovery rate for current xi

    	for (int s=0; s<NoBeta; s++)
		{ 
			if ((s==0) && (xi <= Theta[0]))
                betai=Beta[0];
				
			if ((s==(NoBeta-1)) && (xi>Theta[NoBeta-2]))
                betai=Beta[NoBeta-1];
				
			if ((s>0) && (xi>Theta[s-1]) && (xi<=Theta[s]))
                betai=Beta[s];
		}

		int NoRec=Recovery.size(); 
		int rati;
		for (int l=0; l<NoRec; l++)
		{ 
			if ((l==0) && (xi <= ThresholdR[0]))
			{ 
                recoveryi=Recovery[0];
                rati=rat0;
            } //else	
			if ((l==(NoRec-1)) && (xi>ThresholdR[NoRec-2]))				
			{		 
                recoveryi=Recovery[NoRec-1];
                rati=rat2;
            } //else
			if ((l>0) && (xi>ThresholdR[l-1]) && (xi<=ThresholdR[l]))
			{
                recoveryi=Recovery[l];
                rati=rat1;
			}
            if (recoveryi == 0.6 || rati == rat2) {
                int f = 2;
            }
		}


        // calc conditional prob for all names
		condDefaultProbi=GtoNormalCum((ci-betai*xi-UncondMean)/UncondStdDev);			
		dens=exp(-1*(xi)*(xi)/2)/sqrt(2*3.141592653);
		
		for (int k=0; k<=NoNames; k++)
		{
            // use binormial distribution
            double comb=imsl_d_gamma(NoNames+1)/imsl_d_gamma(k+1)/imsl_d_gamma(NoNames-k+1);
            LossTemp[k]= (1-ProbJ)*dens*comb*pow(condDefaultProbi,k)*pow(1-condDefaultProbi,NoNames-k);
		}

        // dump(LossTemp, "loss temp");
		
		for ( k=0; k<=MaxNoUnits; k++)
		{ 
			if (k % rati ==0 )
			{ 
				CondLoss[k][i] = LossTemp[k/rati];
            } else {
				CondLoss[k][i]=0;
			}		
		}
        // dump(CondLoss, i, "cond loss matrix");
	}

	double wi=1;
	for (int j=0; j<=MaxNoUnits; j++)
	{
		Loss[j] = 0;
        // dump(CondLoss, j, -1, "condloss");
		for (int i=0; i<(NoIntervals+1); i++)
		{ 		
			wi=4/3;

			if ( (i==0) || (i==NoIntervals) ) 
			{wi=27/12;}

			if ( (i==1) || (i==(NoIntervals-1)) ) 
			{wi=0;}

			if ( (i==2) || (i==(NoIntervals-2)) ) 
			{wi=13/12;}

			if (i<=ibound)	
			{
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthL;
			}

			else
			{
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthH;
			}			
		}
    }

    // dump(Loss, "loss");
		
	return Loss;
}

// Calc loss distributions for heterogenous portfolio using RR

// tweak beta
inline double B(double beta, double co)
{
    return beta + co*(1-beta);
}

inline double BLinear(Array<double>& beta, Array<double>& theta, double z)
{
    int j = 1;
    while (j < beta.size() -1 && z > theta[j]) ++j;
    return beta[j-1] + (z-theta[j-1])/(theta[j]-theta[j-1])*(beta[j]-beta[j-1]);
}

Array<double> LossDistRandomNFactorsSimpRR(
	   int					n,
	   Array<double>		&BetaCo, // N beta co
	   Array<double>		&theta, // n-1 theta
       Array<double>		&RecoveryCo,
	   Array<double>		&ThresholdR,
       Array<double>        &Beta,
       Array<double>        &ExpectedRecovery,
	   Array<int>   		&LossGivenDefault,
	   Array<double>		&DefaultProb,
	   double				NoIntervals,
	   double				lbound,							
	   double				ubound,
	   double				mbound,
       int&                 LossUnit,
       bool                 linearBeta)
{
    int NoBeta = BetaCo.size();
    Array<int> LossParam = CalcLossUnit(n, LossGivenDefault, ExpectedRecovery, RecoveryCo);
    int MaxNoUnits = LossParam[1];
    LossUnit = LossParam[0];

    Array<double> Loss(MaxNoUnits+1);

    // calibrate cis for all the names, keep track of average betas 
    Array<double> avgBeta(NoBeta);
    double avgCi(0);

    Array<double> ci(n);
    Array<double> mean(n);
    Array<double> stdDev(n);
    for (int k = 0; k < n; ++k) {

        // calc real betas
        Array<double> beta(NoBeta);
        for (int i = 0; i < NoBeta; ++i) {
            beta[i] = B(Beta[k], BetaCo[i]);
            avgBeta[i] += beta[i];
        }

        // calc unconditional mean and variance
        double UncondMean, UncondStdDev;
   	    UncondMean=beta[0]*GtoNormalDen(theta[0]);
	    UncondStdDev =beta[0]*beta[0]*(GtoNormalCum(theta[0])-theta[0]*GtoNormalDen(theta[0]));	 
	    for (int s=1; s<(NoBeta-1); s++)
	    {
		    UncondMean += beta[s]*(GtoNormalDen(theta[s])-GtoNormalDen(theta[s-1]));
		    UncondStdDev += beta[s]*beta[s]*(GtoNormalCum(theta[s])-GtoNormalCum(theta[s-1])+theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*GtoNormalDen(theta[s]));
	    }
	    UncondMean += beta[NoBeta-1]*(-GtoNormalDen(theta[NoBeta-2]));
	    UncondStdDev += beta[NoBeta-1]*beta[NoBeta-1]*(1-GtoNormalCum(theta[NoBeta-2])+theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2])); 
        UncondStdDev -= UncondMean*UncondMean;
        // fail if variance > 1
	    if (UncondStdDev>1)
        {
            for (int j=0; j<MaxNoUnits+1; j++)
                Loss[j]=-1;
            err("uncond std dev > 1 while calibrating individual names\n");
            return Loss;
        }
  	    UncondStdDev=sqrt(1-UncondStdDev);
        
        mean[k] = UncondMean;
        stdDev[k] = UncondStdDev;

        // bracket ci first;
        GtoNormalCumInv(DefaultProb[k], &ci[k]);
        Array<double> clow(n);
        Array<double> chigh(n);
        Array<double> cmid(n);
        double alpha = 1.05;
        double zero;
        if (ci[k]>0)
        {
            clow[k]=ci[k]/alpha;
            chigh[k]=ci[k]*alpha;
            cmid[k]=(clow[k]+chigh[k])/2;

            // locate upper bound
            double zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
        	while (zeroH<0)
            {
                alpha=1.2*alpha;
                chigh[k]=chigh[k]*alpha;
                zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
            }

            // locate lower bound			
            alpha = 1.05;
        	double zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	
			if (zeroL>0) // what is this? some sort of heuristics?
			{
    	    	alpha=1.2*alpha;
				clow[k]=clow[k]/alpha;
				zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
			}
			if (zeroL>0)
			{  
                clow[k]=-(clow[k]);
                zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 	
				while (zeroL>0)
				{
					alpha=1.2*alpha;
					clow[k]=(clow[k])*alpha;;
					zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 
				}
			}
	  
    		cmid[k]=(clow[k]+chigh[k])/2;
            zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	  
        } else {
            clow[k]=ci[k]*alpha;
            chigh[k]=ci[k]/alpha;
            cmid[k]=(clow[k]+chigh[k])/2;

            double zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  
        	while (zeroL>0)
			{
   				alpha=1.2*alpha;
				clow[k]=clow[k]*alpha;
				zeroL=UncondSurvProbNfact(clow[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k]; 		  			
			}
			
            double zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  	
            alpha=1.5;  
			if (zeroH<0)
			{
    	    	alpha=1.2*alpha;
				chigh[k]=chigh[k]/alpha;
				zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
			}
			if (zeroH<0)
			{  
                chigh[k]=-(chigh[k]);
                zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];	
				while (zeroH<0)
				{
                    alpha=1.2*alpha;
                    chigh[k]=(chigh[k])*alpha;;
    			    zeroH=UncondSurvProbNfact(chigh[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
				}
			}	  
		
		    cmid[k]=(clow[k]+chigh[k])/2;
            zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];  
        }

        // do bisection
		while (fabs(zero)>0.000001)
		{
            if (zero>0)
            {	
                chigh[k]=cmid[k];
                cmid[k]=(chigh[k]+clow[k])/2;
            } else {
                clow[k]=cmid[k];
				cmid[k]=(chigh[k]+clow[k])/2;
			}
            zero=UncondSurvProbNfact(cmid[k],UncondMean,UncondStdDev,beta,theta)-DefaultProb[k];
        }
 
        // found c
 	    ci[k]=cmid[k];
        avgCi += ci[k];
    }

    dump(ci, "ci");

    // calculate r2

    // calc average beta, mean and variance
	for (int i = 0; i < NoBeta; ++i)
        avgBeta[i] /= n;
    avgCi /= n;

    dump(avgBeta, "avgBeta");
    dump(avgCi, "avgCi");

   	double UncondMean=avgBeta[0]*GtoNormalDen(theta[0]);
	double V=avgBeta[0]*avgBeta[0]*(GtoNormalCum(theta[0])-theta[0]*GtoNormalDen(theta[0]));	 
	for (int s=1; s<(NoBeta-1); s++)
	{
		UncondMean=UncondMean+avgBeta[s]*(GtoNormalDen(theta[s])-GtoNormalDen(theta[s-1]));
		V=V+avgBeta[s]*avgBeta[s]*(GtoNormalCum(theta[s])-GtoNormalCum(theta[s-1])+theta[s-1]*GtoNormalDen(theta[s-1])-theta[s]*GtoNormalDen(theta[s]));
	}
	UncondMean=UncondMean+avgBeta[NoBeta-1]*(-GtoNormalDen(theta[NoBeta-2]));
	V=V+avgBeta[NoBeta-1]*avgBeta[NoBeta-1]*(1-GtoNormalCum(theta[NoBeta-2])+theta[NoBeta-2]*GtoNormalDen(theta[NoBeta-2])); 
    V=V-UncondMean*UncondMean;
    // fail if variance > 1
	if (V>1)
    {
        for (int j=0; j<MaxNoUnits+1; j++)
            Loss[j]=-1;
        err("avg variance > 1\n");
        return Loss;
    }
  	double UncondStdDev=sqrt(1-V);        

    int NoRec = RecoveryCo.size();
	double tol=0.01;
	Matrix<double> Param(NoBeta+NoRec-1,3);
	double thresL=ThresholdR[0]; // assume two thresholds
	double thresH=ThresholdR[1];
	double thresM=(thresL+thresH)/2;
	ThresholdR[1]=thresM; // initial guess

        double f  = avgBeta[0];
        f = avgBeta[1];
        f = avgBeta[2];
        f = theta[0];
        f = theta[1];
        f = RecoveryCo[0];
        f = RecoveryCo[1];
        f = RecoveryCo[2];
        f = ThresholdR[0];
        f = ThresholdR[1];

	Param=OrderBetaRecTh(avgBeta,theta,RecoveryCo,ThresholdR);
	double zeroR=UncondSurvProbNBetaNRec(avgCi,UncondMean,UncondStdDev,Param);
	while ((fabs(zeroR-1.0)) >tol && (ThresholdR[1] - ThresholdR[0] > 0.01))
	{
		if (zeroR < 1.0) // expected recovery too small, reduce r2
		{
			thresH=thresM; 
			thresM=(thresL+thresH)/2;	
        } else { // expected recovery too big, increase r2
			thresL=thresM;
			thresM=(thresL+thresH)/2;	
		}	

        ThresholdR[1]=thresM;

        f  = avgBeta[0];
        f = avgBeta[1];
        f = avgBeta[2];
        f = theta[0];
        f = theta[1];
        f = RecoveryCo[0];
        f = RecoveryCo[1];
        f = RecoveryCo[2];
        f = ThresholdR[0];
        f = ThresholdR[1];

        Param=OrderBetaRecTh(avgBeta,theta,RecoveryCo,ThresholdR);
        zeroR=UncondSurvProbNBetaNRec(avgCi,UncondMean,UncondStdDev,Param);
	}
    if (fabs(zeroR-1.0) > tol) { // can not find an appropriate r2
        for (int j=0; j<MaxNoUnits+1; j++)
        {
            Loss[j]=-1;
        }
        if (zeroR > 1.0) {
            err("can not find r2, recovery too big ");
            err(zeroR);
            err("\n");
        }
        else {
            err("can not find r2, recovery too small ");
            err(zeroR);
            err("\n");
        }
        return Loss;
    }

    // dump(ThresholdR, "threshold r");

    // do numerical integration
    double lengthL=2*(mbound-lbound)/NoIntervals;
	double lengthH=2*(ubound-mbound)/NoIntervals;
	int ibound=NoIntervals/2;

    Matrix<double> CondLoss(MaxNoUnits + 1, NoIntervals+1);
    Array<double> LossTemp(MaxNoUnits + 1);
    Array<double> condDefaultProbi(n);
    Array<int> condIndLoss(n);
    
    double xi;
    double betai;
    double reci;
    double dens;
    for (i = 0; i < NoIntervals + 1; ++i) {
    
        //i = 350; // xxx

        // calculate current market factor
		if (i<=ibound)
			xi=lbound+i*lengthL;
		else
			xi=lbound+(NoIntervals/2)*lengthL +lengthH*(i-ibound);
	  
        // calc conditional prob & current loss level for each name under current market factor
        for (int j = 0; j < n; ++j) {
            
            // locate current beta
            if (!linearBeta) {        
                for (int s=0; s<NoBeta; s++)
			    { 
				    if ((s==0) && (xi <= theta[0]))
                        betai=B(Beta[j], BetaCo[0]);
				    
				    else if ((s==(NoBeta-1)) && (xi>theta[NoBeta-2]))
                        betai=B(Beta[j], BetaCo[NoBeta-1]);
				    
				    else if ((s>0) && (xi>theta[s-1]) && (xi<=theta[s]))
                        betai=B(Beta[j], BetaCo[s]);
			    }
            } else {
                betai = BLinear(BetaCo, theta, xi);
                betai = B(Beta[j], betai);
            }

            // calc unconditional prob
            double f1 = mean[j];
            double f2 = stdDev[j];
            double f3 = ci[j];
            double f4 = xi;
            double f5 = betai;
            condDefaultProbi[j] = GtoNormalCum((ci[j] - betai*xi - mean[j])/stdDev[j]);
            double f6 = condDefaultProbi[j];

            // locate current recovery coefficient
            if (xi < ThresholdR[0])
                reci = RecoveryCo[0];
            else if (xi < ThresholdR[1])
                reci = RecoveryCo[1];
            else
                reci = RecoveryCo[2];

            // determine losses
            double indNotional = LossGivenDefault[j]/(1-ExpectedRecovery[j]);
            condIndLoss[j] = indNotional*(1-ExpectedRecovery[j]*reci) + 0.5;
            indNotional = condIndLoss[j];
        }

        // do convolution to calculate loss
        LossTemp = CondLossDist(n, condIndLoss, condDefaultProbi, LossUnit, MaxNoUnits);

        dens = exp(-1*xi*xi/2)/sqrt(2*3.141592653);
        LossTemp *=dens;

        for (int k = 0; k <= MaxNoUnits; ++k)
            CondLoss[k][i] = LossTemp[k];

        // dump(LossTemp, "LossTemp");
        // dump(CondLoss, i, "condloss matrix");
    }

    // do simpson integration
	double wi=1;
	for (int j=0; j<=MaxNoUnits; j++)
	{
		Loss[j] = 0;
        // dump(CondLoss, j, -1, "cond loss");
		for (int i=0; i<(NoIntervals+1); i++)
		{ 		
			wi=4/3;

			if ( (i==0) || (i==NoIntervals) ) 
			{wi=27/12;}

			if ( (i==1) || (i==(NoIntervals-1)) ) 
			{wi=0;}

			if ( (i==2) || (i==(NoIntervals-2)) ) 
			{wi=13/12;}

            // is this right, use 4/3 for all points in the middle?
			if (i<=ibound)
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthL;
			else
				Loss[j] = Loss[j]+wi*CondLoss[j][i]*lengthH;
		}
	}

    // dump(Loss, "Loss");
    return Loss;
}
