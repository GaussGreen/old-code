// Contains the code that derives the implied abse correlation skew
// from the array of expected losses for a set of tranches. In other words allows
// us to express the results of any multiname modeling approach in base correlation
// terms

#include "Magnet/Magnet.h"
#include "Magnet/Function.h"
#include "General/General.h"

#include "Convolution.h"
#include "Quadrature.h"
#include "Tranche.h"

#include "Util.h"

using namespace CM;


// tweaks a scalar beta by some factor 0<c<1 using beta+c(1-beta)
void BetaTweakSingle(double& Beta, double Tweakfact)
{
	Beta += Tweakfact*(1-Beta);
}

// does the same tweak for each element of a vector of Beta's
// Useful for an heterogeenous portfolio
void BetaTweak(Array<double> &Beta, double Tweakfact)
{
	int s=(Beta.size());
	for (int k=0; k<s; k++)
	{
		Beta[k] += Tweakfact*(1-Beta[k]);
	}

}

// Given two vectors of beta tweaks
// the first one by a common factor. The tweaked vector is stored in the second vector
void BetaTweakV(Array<double> &Beta, double	Tweakfact, Array<double> &BetaTweak)
{
	int s=(Beta.size());
	for (int k=0; k<s; k++)
	{
		BetaTweak[k]=Beta[k]+Tweakfact*(1-Beta[k]);
	}

}


// Main Module. Given Expected Losses for a homogeneous beta portfolio
// (same beta, notional, recovery and default probability
// finds the implied base correlation for a given detachment level
// Attachment is a vector  [0,K1,K2,..,1] that partitions total portfolio losses

Array<double> CorrelationSkewSimple(
	 int		   NoNames,
	 int		   NotionalperName, 
	 Array<double> &ExpectedLosses, // Expected Losses in fraction of notional, one per tranche
	 double		   Beta,
	 Array<double> &Attachment,
	 int		   LossperName,
	 double		   DefaultProb)
{ 

	int Notional=NoNames*NotionalperName;
	int NoTranches=ExpectedLosses.size();

	Array<double> Skew(NoTranches-1); // does not include last tranche
	
    int LossUnit=LossperName;
	int MaxNoUnits=NoNames;

   	Array<double> NotTranche(NoTranches);
	Array<double> MarketLoss(NoTranches); // for each tranche [Kl,Ku] denotes the vector of
										  // cumulative expected losses from [0,Ku]
	
    double betasimp=Beta;

	for(int k=0; k<NoTranches; k++)
        MarketLoss[k]=0;
	
	for (k=0; k<(NoTranches-1); k++)
	{	   
        // given current estimated beta, calc loss distribution
        Array<double> lossDist=LossDistConstBetaHomogeneousSpread(NoNames, betasimp, DefaultProb);

        Array<double> TrancheBase(2); // bc tranche
        TrancheBase[0]=0;   
        TrancheBase[1]=Attachment[k+1];

        // given loss distribution calc expected loss for bc tranche
        Array<double> LossModel=ExpectedLossTranche(Notional, LossUnit, TrancheBase, lossDist);

        Array<double> Tranche(2); // actual tranche
        Tranche[0]=Attachment[k];
        Tranche[1]=Attachment[k+1];
        double NotTranche=(Tranche[1]-Tranche[0])*Notional;

        // adjust market loss for all upper tranches
        for (int s=k; s<NoTranches; s++)		
        { 
            MarketLoss[s]=MarketLoss[s]+ExpectedLosses[k]*(Tranche[1]-Tranche[0]);
        }	  

        MarketLoss[k]=MarketLoss[k]/(TrancheBase[1]-TrancheBase[0]);
 	  
        double zeroinitial=0;
        double zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        double count=0; 

        // tweak until Model Expected Losses match Market Expected Losses
        double tweak = 0.4;
        while ((fabs(zeroend)>0.001) && (count <=50))
        { 
            //change sign of tweak when error changes signs
            if ( ((zeroinitial==0)&&(zeroend>0)) || ((zeroinitial*zeroend)<0) )
                tweak=-tweak/3;
		    BetaTweakSingle(betasimp,tweak);
		    lossDist=LossDistConstBetaHomogeneousSpread(NoNames, betasimp, DefaultProb);
		    LossModel=ExpectedLossTranche(Notional, LossUnit,  TrancheBase, lossDist);
		    zeroinitial=zeroend;
		    zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
	
            ++count;
        }
	  
        Skew[k]=betasimp;
	}

	return Skew;
}


// Same as before, but now the portfolio is allowed
// to have different notionals, recoveries and default
// probabilities. Still same beta

Array<double> CorrelationSkew(
	 const int	   Notional, // total notional
	 Array<double> &ExpectedLosses,
	 double		   Beta,
	 Array<double> &Attachment,
	 Array<int>	   &Loss, // how much each name can lose = individual notional * (1-R)
	 Array<double> &DefaultProb)
{ 
	int NoNames=Loss.size();	
	int NoTranches=ExpectedLosses.size();

	Array<double> Skew(NoTranches-1); // does not include last tranche

	Array<int> LossParam=CalcLossUnit(NoNames,Loss);
    int LossUnit=LossParam[0];
	int MaxNoUnits=LossParam[1];

	Array<double> NotTranche(NoTranches);
	Array<double> MarketLoss(NoTranches); // for each tranche [Kl,Ku] denotes the vector of
										  // cumulative expected losses from [0,Ku]
	
	for(int k=0; k<NoTranches; k++)
        MarketLoss[k]=0;	

    double betasimp=Beta;    
	for (k=0; k<(NoTranches-1); k++)
	{
        double tweak=0.4;
      
        // calc loss distribution
        Array<double> lossDist=LossDistConstBeta(NoNames, betasimp, Loss, DefaultProb);

        Array<double> Tranche(2);    // Actual tranche traded in the market [Kl,Ku]
        Tranche[0]=Attachment[k];
        Tranche[1]=Attachment[k+1];  

        Array<double> TrancheBase(2); // Tranche that defines BC at dettachment Ku, [0,Ku]    
        TrancheBase[0]=0;	
        TrancheBase[1]=Attachment[k+1];

        double NotTranche=(Tranche[1]-Tranche[0])*Notional;
	  
        for (int s=k; s<NoTranches; s++)
            MarketLoss[s]=MarketLoss[s]+ExpectedLosses[k]*(Tranche[1]-Tranche[0]);  

        MarketLoss[k]=MarketLoss[k]/(TrancheBase[1]-TrancheBase[0]);
	  
        Array<double> LossModel=ExpectedLossTranche(Notional, LossUnit, TrancheBase, lossDist);
 	  
        double zeroinitial=0;      
        double zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        double count=0; 

        while ((fabs(zeroend)>0.001) && (count <=50))
        { 
    		//change sign of tweak when error changes signs
		    if ( ((zeroinitial==0)&&(zeroend>0)) || ((zeroinitial*zeroend)<0) )
    		    tweak=-tweak/3;
            BetaTweakSingle(betasimp,tweak);
            lossDist=LossDistConstBeta(NoNames, betasimp, Loss, DefaultProb);
            LossModel=ExpectedLossTranche(Notional, LossUnit,  TrancheBase, lossDist);
            zeroinitial=zeroend;
            zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
	
            count =count+1;
        }
	  
        Skew[k]=betasimp;
	}

	return Skew;
}


// Perform the same routine as before with an heterogenous
// portfolio. The same tweak factor is applied to each names's beta

Array<double> CorrelationSkewGeneral(
	 const int	   Notional,
	 Array<double> &ExpectedLosses,
	 Array<double> &Beta,
	 Array<double> &Attachment,
	 Array<int>	   &Loss,
	 Array<double> &DefaultProb)
{ 
	int NoNames=Beta.size();
	int NoTranches=ExpectedLosses.size();

	Array<double> Skew(NoTranches-1);

	Array<int> LossParam=CalcLossUnit(NoNames,Loss);
    int LossUnit=LossParam[0];
	int MaxNoUnits=LossParam[1];

	Array<double> MarketLoss(NoTranches);
	Array<double> NotTranche(NoTranches);    

	for(int k=0; k<NoTranches; k++)
        MarketLoss[k]=0;	
	
	for (k=0; k<(NoTranches-1); k++)
	{
        double tweak=0.4;
      
        Array<double> lossDist=LossDist(NoNames, Beta, Loss, DefaultProb);

        Array<double> Tranche(2);     
        Tranche[0]=Attachment[k];	  
        Tranche[1]=Attachment[k+1];  

        Array<double> TrancheBase(2);
        TrancheBase[0]=0;
        TrancheBase[1]=Attachment[k+1];

        double NotTranche=(Tranche[1]-Tranche[0])*Notional;
	  
        for (int s=k; s<NoTranches; s++)
        {
            MarketLoss[s]=MarketLoss[s]+ExpectedLosses[k]*(Tranche[1]-Tranche[0]);
        }	  

        MarketLoss[k]=MarketLoss[k]/(TrancheBase[1]-TrancheBase[0]);
	  
        Array<double> LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
 	  
        double zeroinitial=0;
        double zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        double count=0; 

        while ((fabs(zeroend)>0.001) && (count <=50))
        { 
	        if ( ((zeroinitial==0)&&(zeroend>0)) || ((zeroinitial*zeroend)<0) )
                tweak=-tweak/3;

            BetaTweak(Beta,tweak);
            lossDist=LossDist(NoNames, Beta, Loss, DefaultProb);
            LossModel=ExpectedLossTranche(Notional, LossUnit,TrancheBase, lossDist);
            zeroinitial=zeroend;
            zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
	
            count =count+1;
        }

        double temp=0;
        for (int i=0; i<NoNames; i++)
        {
            temp=temp+Beta[i];
        } 
	  
        Skew[k]=temp/NoNames; // average beta?
	}

	return Skew;
}



// Same as before, but now the portfolio is allowed
// to have different notionals, recoveries and default
// probabilities. Still same beta

Array<double> CorrelationSkew2(
	 const int	   Notional, // total notional
	 Array<double> &ExpectedLosses,
	 double		   Beta,
	 Matrix<double> &Attachment,
	 Array<int>	   &Loss, // how much each name can lose = individual notional * (1-R)
	 Array<double> &DefaultProb)
{ 
	int NoNames=Loss.size();	
	int NoTranches=ExpectedLosses.size();

	Array<double> Skew(NoTranches);

	Array<int> LossParam=CalcLossUnit(NoNames,Loss);
    int LossUnit=LossParam[0];
	int MaxNoUnits=LossParam[1];

	Array<double> NotTranche(NoTranches);
	Array<double> MarketLoss(NoTranches);
	
    err("start calc correlation skew\n");

	for (int k=0; k < NoTranches; k++)
	{
        err("calc tranche # "); err(k); err("\n");
        double tweak=0.4;
      
        // calc loss distribution
        Array<double> lossDist=LossDistConstBeta(NoNames, Beta, Loss, DefaultProb);

        Array<double> Tranche(2);    // Actual tranche traded in the market [Kl,Ku]
        Tranche[0]=Attachment[k][0];
        Tranche[1]=Attachment[k][1];  
	  
        Array<double> LossModel=ExpectedLossTranche(Notional, LossUnit, Tranche, lossDist);
 	  
        double difference = ExpectedLosses[k]-LossModel[0];

        double count=0; 

        double bHi = 1;
        double bLow = 0;
        double bMid = Beta;

        while ((fabs(difference)>0.001) && (count <=50))
        {
            if (difference > 0) {
                // model loss too small, reduce beta
                bHi = bMid;
                bMid = (bLow + bMid)/2;
            } else {
                // model loss too big, increase beta
                bLow = bMid;
                bMid = (bLow + bHi)/2;
            }

            lossDist = LossDistConstBeta(NoNames, bMid, Loss, DefaultProb);
            LossModel = ExpectedLossTranche(Notional, LossUnit,  Tranche, lossDist);

            difference = ExpectedLosses[k]-LossModel[0];
            count = count+1;
        }
        Skew[k] = bMid;
	}

    err("end calc skew\n");

	return Skew;
}


// Perform the same routine as before with an heterogenous
// portfolio. The same tweak factor is applied to each names's beta

Array<double> CorrelationSkewGeneral2(
	 const int	   Notional,
	 Array<double> &ExpectedLosses,
	 Array<double> &Beta,
	 Matrix<double> &Attachment,
	 Array<int>	   &Loss,
	 Array<double> &DefaultProb)
{ 
	int NoNames=Beta.size();
	int NoTranches=ExpectedLosses.size();

	Array<double> Skew(NoTranches);

	Array<int> LossParam=CalcLossUnit(NoNames,Loss);
    int LossUnit=LossParam[0];
	int MaxNoUnits=LossParam[1];

	Array<double> MarketLoss(NoTranches);
	Array<double> NotTranche(NoTranches);    

	for(int k=0; k < NoTranches; k++)
        MarketLoss[k]=0;	
	
	for (k=0; k < NoTranches; k++)
	{
        double tweak=0.4;
      
        Array<double> lossDist=LossDist(NoNames, Beta, Loss, DefaultProb);

        Array<double> Tranche(2);     
        Tranche[0]=Attachment[k][0];	  
        Tranche[1]=Attachment[k][1];  
	  
        Array<double> LossModel=ExpectedLossTranche( Notional, LossUnit, Tranche, lossDist);
 	  
        double zeroinitial=0;
        double zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        double count=0; 

        while ((fabs(zeroend)>0.001) && (count <=50))
        { 
	        if ( ((zeroinitial==0)&&(zeroend>0)) || ((zeroinitial*zeroend)<0) )
                tweak=-tweak/3;

            BetaTweak(Beta,tweak);
            lossDist=LossDist(NoNames, Beta, Loss, DefaultProb);
            LossModel=ExpectedLossTranche(Notional, LossUnit,Tranche, lossDist);
            zeroinitial=zeroend;
            zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
	
            count =count+1;
        }

        double temp=0;
        for (int i=0; i<NoNames; i++)
        {
            temp=temp+Beta[i];
        } 
	  
        Skew[k]=temp/NoNames; // average beta?
	}

	return Skew;
}




// Calculates Correlation Skew for a General Heterogenous
// Portfolio using a Secant False Position Method
// See Numerical Recipes in C for Details 

Array<double> CorrelationSkewGeneralSecant(
	 const int	   Notional,
	 Array<double> &ExpectedLosses,
	 Array<double> &Beta,
	 Array<double> &Attachment,
	 Array<int>	   &Loss,
	 Array<double> &DefaultProb)
{ 
	int NoNames=Beta.size();	
	int NoTranches=ExpectedLosses.size();

	Array<double> Skew(NoTranches-1);

	Array<double> BetaL(NoNames);
	Array<double> BetaH(NoNames);
	Array<double> BetaTemp(NoNames);
	
	Array<int> LossParam=CalcLossUnit(NoNames,Loss);
    int LossUnit=LossParam[0];
	int MaxNoUnits=LossParam[1];

	Array<double> MarketLoss(NoTranches);
	Array<double> NotTranche(NoTranches);    

	for(int k=0; k<NoTranches; k++)
        MarketLoss[k]=0;
	
	
	for (k=0; k<(NoTranches-1); k++)
	{	   
        Array<double> BetaTemp(NoNames);

        // looks for a lower bound to bracket Base Correlation for Tranche K
        // keep reducing beta, thus reducing modelLoss

    	double tweakL=-0.4;
      
        Array<double> lossDist; // =LossDist(NoNames, Beta, Loss, DefaultProb); not necessary

        Array<double> Tranche(2);
        Tranche[0]=Attachment[k];
        Tranche[1]=Attachment[k+1];

        Array<double> TrancheBase(2);        	  
        TrancheBase[0]=0;	
        TrancheBase[1]=Attachment[k+1];

        double NotTranche=(Tranche[1]-Tranche[0])*Notional;
	  
        for (int s=k; s<NoTranches; s++)
        { 
            MarketLoss[s]=MarketLoss[s]+ExpectedLosses[k]*(Tranche[1]-Tranche[0]);
        }
	  
        MarketLoss[k]=MarketLoss[k]/(TrancheBase[1]-TrancheBase[0]);
	  
        BetaTweakV(Beta,tweakL,BetaL);
        lossDist=LossDist(NoNames, BetaL, Loss, DefaultProb);

        Array<double> LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
 	      
        double zeroendL=(MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        while (zeroendL<0)
        { 
            tweakL=1.1*tweakL;
            BetaTweakV(Beta,tweakL,BetaL);
            lossDist=LossDist(NoNames, BetaL, Loss, DefaultProb);
            LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
            zeroendL= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]); 
        }
	
        // looks for an upper bound to bracket Base Correlation for Tranche K
        // keep increasing beta, thus increasing modelloss
        double tweakH=0.4;
 
        BetaTweakV(Beta,tweakH,BetaH);
        lossDist=LossDist(NoNames, BetaH, Loss, DefaultProb);
	  
        LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
	 
        double zeroendH= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
	 
        while (zeroendH>0)
        { 
            tweakH=1.1*tweakH;
            BetaTweakV(Beta,tweakH,BetaH);
            lossDist=LossDist(NoNames, BetaH, Loss, DefaultProb);
            LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
            zeroendH= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        }

        // uses secant method to find next trial point updated lower and upper brackets
        // so that the zero of the equation stays within the brackets
	  
        // find the line intersection
        double tweak=tweakL-zeroendL*(tweakH-tweakL)/(zeroendH-zeroendL);
        BetaTweakV(Beta,tweak,BetaTemp);
        lossDist=LossDist(NoNames, BetaTemp, Loss, DefaultProb);
        LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
	
        double zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
		double count=0;

        while ((fabs(zeroend)>0.001) && (count <=50))
        { 
            count=count+1;

            if ( zeroend>0 ) // ? should this be > 0, or < 0
            { 
                zeroendH=zeroend;
                BetaH=BetaTemp;
                tweakH=tweak;
                tweak=tweakL-zeroendL*(tweakH-tweakL)/(zeroendH-zeroendL);
            } else {
                zeroendL=zeroend;
                BetaL=BetaTemp;
                tweakL=tweak;
                tweak=tweakL-zeroendL*(tweakH-tweakL)/(zeroendH-zeroendL);
            }
		
            BetaTweakV(Beta,tweak,BetaTemp);
            lossDist=LossDist(NoNames, BetaTemp, Loss, DefaultProb);
            LossModel=ExpectedLossTranche( Notional, LossUnit, TrancheBase, lossDist);
            zeroend= (MarketLoss[k]-LossModel[0])/(MarketLoss[k]);
        }  
	  
        // return Base Correlation
        double temp=0;
        for (int i=0; i<NoNames; i++)
        {
            temp=temp+BetaTemp[i];
        } 
	  
        Skew[k]=temp/NoNames;
	}

	return Skew;
}
