
// This module contains code that calculates Expected Losses and Tranche
// Spreads given information about the unconditional loss distribution
// for an hetereogeneous portfolio of names


#include "Magnet/Magnet.h"
#include "General/General.h"

#include "Convolution.h"
#include "Quadrature.h"
#include "Tranche.h"

#include "Util.h"

using namespace CM;


// Simple routine that calculates the losses on a tranche given
// the attachment and detachment points as well as the actual portfolio
// losses following EL=max(0,min(x-Kl,Ku-Kl))

double Collar(
	 const double Loss,
	 const double TrancheLow,
	 const double TrancheHigh)

{
	if (Loss < TrancheLow)
		return 0;

    if (Loss < TrancheHigh)
        return Loss - TrancheLow;

    return TrancheHigh-TrancheLow;
}


// Calculates the percentual expected losses for a complete set
// of tranches, given the total notional of the portfolio, the
// array of attachment points and the unconditional portfolio loss
// distribution

// Note: the attachment vector includes all attachment points with 0 being
// the initial component and 1 its final component

Array<double> ExpectedLossTranche(
	 const double  Not,
	 const int     LossUnit,
	 Array<double> &Attachment, 
	 Array<double> &LossDist)
{	
	int MaxNoUnits=LossDist.size();
	int NoTranches=(Attachment.size())-1;

	// Calculate Tranche Prices
	Array<double> CumExpectedLoss(NoTranches);
	if (LossDist[0]==-1)
	{ 
		for (int j=0; j<NoTranches; ++j)
            CumExpectedLoss[j]=1;

		return CumExpectedLoss;	
	}
	
	// intializes expected losses
	for (int k=0; k<NoTranches; ++k)
	{ 
		CumExpectedLoss[k]=0;
	}

	// Loops over all loss probabilities of the distribution	
	for (int j=0; j<=MaxNoUnits; j++)
	{
		// for each element of the distribution updates all tranches
		for(int k=0; k<NoTranches; k++)
		{
            CumExpectedLoss[k] += 
                Collar(j*LossUnit,Attachment[k]*Not,Attachment[k+1]*Not)*LossDist[j];				
		}
	}
    // dump(CumExpectedLoss, "expected loss");

	// divides by total notional of the tranche
	for( k=0; k<NoTranches; k++)
	{ 
        CumExpectedLoss[k]=CumExpectedLoss[k]/(Not*(Attachment[k+1]-Attachment[k]));		
	}

    // dump(CumExpectedLoss, "expected loss");
	return CumExpectedLoss;
} 

Array<double> ExpectedLossTranchelet(
	 const double	  Not,
	 const int		  LossUnit,
	 Matrix<double>	  &Attachment,
	 Array<double>	  &LossDist)
{
    int col = Attachment.colSize();
    int NoTranches = Attachment.rowSize();

    Array<double> CumExpectedLoss(NoTranches);
    for (int i = 0; i < NoTranches; ++i)
        CumExpectedLoss[i] = 0;

    if (col != 2) {
        err("matrix of attachments must have col size 2");
        return CumExpectedLoss;
    }

  	int MaxNoUnits=LossDist.size();

	// Loops over all loss probabilities of the distribution	
	for (int j=0; j<=MaxNoUnits; j++)
	{
		// for each element of the distribution updates all tranches
		for(int k=0; k<NoTranches; k++)
		{
            CumExpectedLoss[k] += 
                Collar(j*LossUnit,Attachment[k][0]*Not,Attachment[k][1]*Not)*LossDist[j];				
		}
	}
    // dump(CumExpectedLoss, "expected loss");

	// divides by total notional of the tranche
	for(int k=0; k<NoTranches; k++)
	{ 
        CumExpectedLoss[k]=CumExpectedLoss[k]/(Not*(Attachment[k][1]-Attachment[k][0]));		
	}

    // dump(CumExpectedLoss, "expected loss");

    return CumExpectedLoss;
}

Array<double> ExpectedLossTrancheBuck(
	 double		   Not,
	 double		   LossUnit,
	 int			MaxNoUnits,
	 Array<double> &Attachment, 
	 Matrix<double> &LossDist)
{
	int NoTranches=(Attachment.size())-1;

	// Calculate Tranche Prices	
	Array<double> CumExpectedLoss (NoTranches);		
		
	if (LossDist[0][0]==-1)
	{ 
		for (int j=0; j<NoTranches; j++)
            CumExpectedLoss[j]=1;
		return CumExpectedLoss;	
	}

	// intializes expected losses
	for (int k=0; k<NoTranches; k++)
	{ 
		CumExpectedLoss[k]=0;
	}
	
	// Loops over all loss probabilities of the distribution
	for (int j=0; j<=MaxNoUnits; j++)
	{		
		// for each element of the distribution updates all tranches
		for(int k=0; k<NoTranches; k++)
            CumExpectedLoss[k] += 
                Collar(LossDist[j][1]*LossUnit,Attachment[k]*Not,Attachment[k+1]*Not)*LossDist[j][0];		
	}

	// divides by total notional of the tranche
	for( k=0; k<NoTranches; k++)
	{ 
      CumExpectedLoss[k]=CumExpectedLoss[k]/(Not*(Attachment[k+1]-Attachment[k]));		
	}
	
	return CumExpectedLoss;
} 


// Does the same for a multiperiod model. 
// Interpolation needs to be improved. Do not use yet.

Matrix<double> SpreadTranche(
	 const double  Not,
	 const int     LossUnit,
	 const int	   MaxNoUnits,
	 Array<double> &Attachment,
	 Matrix<double> &LossDist,
	 Array<double>	&Maturity,
	 double			time)
{
	int NoTranches=(Attachment.size())-1;
	
	int NoPeriods=LossDist.colSize();

	// Calculate Tranche Prices
	Matrix<double> CumExpectedLoss(NoTranches,NoPeriods);		
	Matrix<double> TrancheSpread(NoTranches,NoPeriods);		

	if (LossDist[0][0]==-1)
	{ 
		for (int l=0; l<NoPeriods; l++)
		{
		    for (int j=0; j<NoTranches; j++)
    		{
                CumExpectedLoss[j][l]=1;
            }
		}
		return CumExpectedLoss; // XXX, should return dummy spread matrix, not expected loss
	}

	// initializes expected losses
    for (int k=0; k<NoTranches; ++k)
    {
        for (int l=0; l<NoPeriods; ++l)
        {
            CumExpectedLoss[k][l]=0;
        }
    }

    // iterate over the periods
	for (int l=0; l<NoPeriods; l++)
	{
        // calc expected tranche loss for this period
		for (int j=0; j<=MaxNoUnits; j++)
		{
			for(int k=0; k<NoTranches; k++)
			{
                CumExpectedLoss[k][l]+=
                    Collar(j*LossUnit,Attachment[k]*Not,Attachment[k+1]*Not)*LossDist[j][l]
                        /(double(Not)*(Attachment[k+1]-Attachment[k]));						
			}
		}

    	double num; double denom;

        // calc spread for each tranche
	    for (k = 0; k < NoTranches; k++)
		{ 
			num=CumExpectedLoss[k][l];
			denom=0;
	
			double spreadi;

            // look at all previous periods and calc denominators
            // this logic seem to be wrong
			for (int s=0; s<=l; s++)
			{
				int n;
				
                // calc raw spread
				if (s==0)			
				{
                    n=Maturity[s]/time; // day count
                    spreadi = (1-CumExpectedLoss[k][s]);
                    spreadi = -log(spreadi)/(Maturity[s]);
                } else {
                    n=(Maturity[s]-Maturity[s-1])/time; // day count
                    spreadi = (1-CumExpectedLoss[k][s])/(1-CumExpectedLoss[k][s-1]);
                    spreadi = -log(spreadi)/(Maturity[s]-Maturity[s-1]);
				}
						
				for (int y=0; y<n; y++)
				{
					if (s==0)
					{
                        denom=denom+exp(-spreadi*(y+1)*time);
                    } else {
                        denom=denom+(1-CumExpectedLoss[k][s-1])*exp(-spreadi*(y+1)*time);
					}
				}
			}

			if (k==0)
			{
                // for equity tranche
				TrancheSpread[k][l]= (num-denom*500/10000)/(Attachment[1]-Attachment[0]);
            } else {
				TrancheSpread[k][l]=(num/denom)*10000;
			}
		}
	}
	
	return TrancheSpread;
}






