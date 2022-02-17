
// This module contains the code that allows to find the loss distribution of
// an heterogeneous portfolio of different names given that the default probabilities
// are independet



#include "Magnet/Magnet.h"
#include "General/General.h"

#include "Convolution.h"
#include "Util.h"

using namespace CM;


//calculates the Greatest Common Divisor of two given numbers
int GCD(int x,int y)
{
    while( x != 0 ) {
        int oldY = y;
        y = x;
        x = oldY % x;
    }
    return y;
}


//calculates the Least Common Multiple of two given numbers
int LCM(int x,int y)
{
    return x*y/(GCD(x,y));
}


// finds the loss unit and max no of units from a given array of loss given defaults
// Loss Unit is taken as the Greatest Common Divisor of all possible name losses
// Max No Units is the maximum portfolio loss expressed in loss units
Array<int> CalcLossUnit(int	NoNames, const Array<int> LossGivenDefault)
{
	Array<int> LossUnit(2);
	
	LossUnit[0]=LossGivenDefault[0];
   
    int maxUnit=0; 
 
    // calculates the GCD of all names and adds total no. of units
    for(int i=0; i<NoNames; i++)
	{
        maxUnit= maxUnit+LossGivenDefault[i];
        if (i<(NoNames-1))
        { 
		  LossUnit[0]=GCD(LossUnit[0],LossGivenDefault[i+1]);
        }
	}
	
	maxUnit/=LossUnit[0];
	LossUnit[1]=maxUnit;
	return LossUnit;
}

Array<int> CalcLossUnit(int	NoNames, 
                        const Array<int>& LossGivenDefault,
                        const Array<double>& ExpectedRecovery,
                        const Array<double>& Co)
{
    int size = LossGivenDefault.size();
    Array<int> Lgd(size);
    Array<int> Lou(2);
    Lou[0] = LossGivenDefault[0];
    int maxunit;
    for (int i = Co.size()-1; i >= 0; --i) {
        maxunit = 0;
        for (int j = 0; j < size; ++j) {
            double indNotional = LossGivenDefault[j]/(1-ExpectedRecovery[j]);
            Lgd[j] = indNotional*(1-ExpectedRecovery[j]*Co[i]) + 0.5;
            maxunit += Lgd[j];
        }
        Array<int> current = CalcLossUnit(NoNames, Lgd);
        if (i == Co.size() - 1)
            Lou[0] = current[0];
        else
            Lou[0] = GCD(Lou[0], current[0]);
    }
    Lou[1] = maxunit/Lou[0];

    int f = Lou[0];
    int g = Lou[1];

	return Lou;
}


// Calculates the loss Distribution of the portfolio conditioned
// on the realization of the market factor by implementing a Standard
// convolution routine

Array<double> CondLossDist(
	 int					NoNames,
	 const Array<int>       &LossGivenDefault,
	 const Array<double>	&DefaultProb,
     int                    LossUnit,
     int                    maxNoUnits)
{
    //dump(LossGivenDefault, "LossGovenDefault");
    //dump(DefaultProb, "DefaultProb");

    // some preparation work
    Array<int> UnitLossVect(NoNames); // number of loss units given default for name i
	Array<int> maxLossUnits(NoNames); // for each name i defines the sum of loss units for the previous i-1 names
	double temp; double cum=0;
	maxLossUnits[0]=0;
	for(int k=0; k<NoNames; k++)
	{
        temp = (LossGivenDefault[k])/LossUnit;
        UnitLossVect[k]=temp;
        cum=cum+temp;
        if (k<(NoNames-1))
		{ 
		    maxLossUnits[k+1]=cum; 
		}
	}
	  	     
  
    // initialization
	Array<double> PDFLoss(maxNoUnits+1);
	Array<double> PDFLosstemp(maxNoUnits+1);
 
	PDFLoss[0]=1;
	PDFLosstemp[0]=1;
	
	for(int s=1; s<= maxNoUnits; s++)	
	{
		PDFLoss[s]=0;
		PDFLosstemp[s]=0;
	}
      
    // Loop over names
    for (int i=0; i<NoNames; i++)
	{
		for (int j=0; j<= maxLossUnits[i]; j++)
        // Loop over all possible units and add survival probabilities
		{
			int uj = j + UnitLossVect[i];
			PDFLosstemp[j]=PDFLosstemp[j]-PDFLoss[j]*(DefaultProb[i]); //update survival probability
			PDFLosstemp[uj]=PDFLosstemp[uj]+PDFLoss[j]*(DefaultProb[i]); // update default probability
		}     
        PDFLoss=PDFLosstemp;
    } 

    //dump(PDFLoss, "PDFLoss");
    
    return PDFLoss;
}



Array<double> CondLossDist(
	 int					NoNames,
	 const Array<int>       &LossGivenDefault,
	 const Array<double>	&DefaultProb)
{
	Array <int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
    return CondLossDist(NoNames, LossGivenDefault, DefaultProb, LossParam[0], LossParam[1]);
}

// Finds the conditional loss distribution for a given array of maturity dates
// given a term structure of default probalities for each date and loss given
// default for each name. The code is almost identical to the previous one but,
// instead, the output is a matrix with each colummn corresponding to a loss distribution
// at a fixed maturity date

Matrix<double> CondLossDistMultiper(
	 const int				NoNames,
	 const Array<int>		&LossGivenDefault,
	 const Matrix<double>	&DefaultProb)
{
	int NoPeriods=DefaultProb.colSize();

	Array <int> LossParam=CalcLossUnit(NoNames,LossGivenDefault);
    int LossUnit=LossParam[0];
	int maxNoUnits=LossParam[1];
  
    // some preparation work
    Array<int> UnitLossVect(NoNames); 
	Array<int> maxLossUnits(NoNames);
	maxLossUnits[0]=0;
	double temp; double cum=0;	
	for(int k=0; k<NoNames; k++)
	{
		temp = (LossGivenDefault[k])/LossUnit;	
        UnitLossVect[k] = temp;
		cum=cum+temp;
        if (k<(NoNames-1))
		{
            maxLossUnits[k+1]=cum;
		}
	}
	  	     
    // initialization
 	Matrix<double> PDFLoss(maxNoUnits+1,NoPeriods);
	Array<double> PDFLosstemp(maxNoUnits+1, NoPeriods);	
	for (int m=0; m<NoPeriods; m++)
	{	
        PDFLoss[0][m]=1;
	}

 	for(int s=1; s<= maxNoUnits; s++)	
	{
        for (m=0; m<NoPeriods; m++)
        {
            PDFLoss[s][m]=0;
        }	
	}
   
    // Loop over names
	for (m=0; m<NoPeriods; m++)
	{
		
        PDFLosstemp[0]=1;
        for (s=1; s<=maxNoUnits; s++)
        {
            PDFLosstemp[s]=0;
        } 
			
        for(int i=0; i<NoNames; i++)
		{
			// Loop over all possible units and add survival probabilities
	   
			for (int j=0; j<= maxLossUnits[i]; j++)
			{
				int uj = j + UnitLossVect[i];
				PDFLosstemp[j]=PDFLosstemp[j]-PDFLoss[j][m]*(DefaultProb[i][m]);
				PDFLosstemp[uj]=PDFLosstemp[uj]+PDFLoss[j][m]*(DefaultProb[i][m]);
			}
     
			for (int l=0; l<NoNames; l++)
			{
			  PDFLoss[l][m]=PDFLosstemp[l];
			}
        }
    }

    return PDFLoss;
}


// calculate loss distribution using Hull-White probability bucketing
// this function does not seem to do the right thing, do not use

Matrix<double> CondLossDistUnits(
	 int					TotalLossUnits,
	 Array<double>			&LossUnitsGivenDefault,
	 Array<double>			&DefaultProb)
{

	int NoNames=DefaultProb.size();

    // some preparation work
	Array<double> maxLossUnits(NoNames); // for each name i defines the sum of loss units for the previous i-1 names
	double temp; double cum=0;
	maxLossUnits[0]=0;
	for(int k=0; k<NoNames; k++)
	{
        temp = (LossUnitsGivenDefault[k]);
        cum=cum+temp;
        if (k<(NoNames-1))
        { 
            maxLossUnits[k+1]=cum; 
		}
	}     
	double maxNoUnits=cum;

    // initialization
	Array<double> PDFLoss(TotalLossUnits+1);
	Array<double> PDFLosstemp(TotalLossUnits+1);
 
	Array<double> Mean(TotalLossUnits+1);
	Array<double> MeanTemp(TotalLossUnits+1);

    // with 0 names, no loss prob 1, mean loss 0
	Mean[0]=0;
	MeanTemp[0]=0;
	PDFLoss[0]=1;
	PDFLosstemp[0]=1;
	
	for(int s=1; s<= TotalLossUnits; s++)
	{
		PDFLoss[s]=0;
		PDFLosstemp[s]=0;
		Mean[s]=0.5*(2*s-1);		
		MeanTemp[s]=0.5*(2*s-1);
	}
     
   	int uj;
    double tol=0.00000000001;

    // loop over names
    for (int i=0; i<NoNames; i++)
	{    	
	    // loop over all units and add survival probabilities
		for (int j=0; j<= ceil(maxLossUnits[i]); j++)  
		{
			if (j==0)
			{ 
				uj=floor(LossUnitsGivenDefault[i]-tol)+1;	
            } else {
                double jump=MeanTemp[j]+LossUnitsGivenDefault[i];
				uj=floor(jump-tol)+1;
			}
			
			PDFLosstemp[j]=PDFLosstemp[j]-PDFLoss[j]*(DefaultProb[i]); //update survival probability
			double c=PDFLosstemp[uj]+PDFLoss[j]*(DefaultProb[i]); // update default probability

            // is this condition right?
			if (c>0)
			{
				MeanTemp[uj]=
                    (MeanTemp[uj]*PDFLosstemp[uj] + 
                        (MeanTemp[j]+LossUnitsGivenDefault[i])*PDFLoss[j]*DefaultProb[i])/c;
            } else {
				MeanTemp[uj]=MeanTemp[uj];
			}
			
		    PDFLosstemp[uj]=c;		
		}

        PDFLoss=PDFLosstemp;
        Mean=MeanTemp;
	} 
	
  	Matrix<double> Dist(TotalLossUnits+1,2);
	for (int h=0; h<=TotalLossUnits; h++)
	{ 
		Dist[h][0]=PDFLoss[h];
		Dist[h][1]=Mean[h];
	}
		
	return Dist; 
}



