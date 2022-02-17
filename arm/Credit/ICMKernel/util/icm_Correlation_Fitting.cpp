#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Correlation_Fitting.CPP
	PROJECT:	UTIL
	
	DESCRIPTION:	Containor for a Sector Correlation

  -----------------------------------------------------------------

 	CREATION:	February, 2006

	LAST MODIF:	February, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */


#include "ICMKernel\util\icm_correlation_fitting.h"

# include "ICMKernel/util/icm_BFGS.h"
# include "ICMKernel/glob/icm_maths.h"


void	ICM_Correlation_Fitting::Compute_BestFit_CorrelationMatrix()
{
	// ---------------------------------------------------------------------
	// TEST
	if (its_CorrelationMatrix == NULL)
		ICMTHROW(ERR_INVALID_DATA,"Null matrix in ICM_Correlation_Fitting!");
	// ---------------------------------------------------------------------

	its_CorrelationSize	=	its_CorrelationMatrix->GetSize();
	// ---------------------------------------------------------------------
	// SOLVER Parameters

	bool	res;
	double	val;
	int		nbiter;

	int		TheDim;
	int		i;

	// ---------------------------------------------------------------------
	double*	TheX	=	NULL;

	its_Steps.clear();
	its_NbMaxIter = 100;
	nbiter	=	0;
	val		=	0.0;

	its_correl_parameters.clear();
	its_Steps.clear();

	string	Correl_Type_String;

	switch (its_Correlation_Fit_Type)
	{
		case	qCORREL_1F_SINGLE_VALUE:
			
			Correl_Type_String	=	"ONE FACTOR SINGLE VALUE";
			
			its_dim	=	1;
			its_1F_Single_Value = its_CorrelationMatrix->Average();

			break;

		case	qCORREL_1F_BETAS:
			
			Correl_Type_String	=	"ONE FACTOR BETAS";
			
			// Dimension
			TheDim	=	its_CorrelationMatrix->GetSize();
			its_dim	=	TheDim;

			its_correl_parameters.resize(its_dim);
			its_Steps.resize(its_dim);

			if (TheX)
				delete	TheX;
			TheX	=	new double [TheDim];

			for (i=0; i<its_dim; i++)
			{
				TheX[i]	=	0.5;
				its_Steps[i]		=	0.01;
			}

			res = BFGS(	
					ff1::mem_call(&ICM_Correlation_Fitting::TheBetaBestFitFunction, *((ICM_Correlation_Fitting*)this)) , 
					ff1::mem_call(&ICM_Correlation_Fitting::TheBetaBestFitGradientFunction, *((ICM_Correlation_Fitting*)this))
					).Minimize(TheX, TheDim, 0.000001, &nbiter, &val );	

			its_1F_Betas.clear();
			its_1F_Betas.resize(TheDim);

			for (i=0; i<TheDim; i++)
				its_1F_Betas[i]	=	its_correl_parameters[i];

			break;


		case	qCORREL_2F_SAME_INTER_SAME_INTRA:

			Correl_Type_String	=	"TWO FACTORS SAME INTER SAME INTRA";

			// Dimension
			TheDim	=	2;
			its_dim	=	TheDim;

			its_correl_parameters.resize(its_dim);
			its_Steps.resize(its_dim);

			if (TheX)
				delete	TheX;
			TheX	=	new double [TheDim];

			for (i=0; i<its_dim; i++)
			{
				TheX[i]	=	0.5;
				its_Steps[i]		=	0.01;
			}

			res = BFGS(	
					ff1::mem_call(&ICM_Correlation_Fitting::The2FSameInterSameIntraBestFitFunction, *((ICM_Correlation_Fitting*)this)) , 
					ff1::mem_call(&ICM_Correlation_Fitting::The2FSameInterSameIntraBestFitGradientFunction, *((ICM_Correlation_Fitting*)this))
					).Minimize(TheX, TheDim, 0.000001, &nbiter, &val);	

			its_2F_Inter_Sector_Value	=	its_correl_parameters[0];
			its_2F_Intra_Sector_Value	=	its_correl_parameters[1];

			break;


		case	qCORREL_2F_SAME_INTER_DIFF_INTRA:
			Correl_Type_String	=	"TWO FACTORS SAME INTER DIFF INTRA";

			ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
			break;

		case	qCORREL_2F_DIFF_INTER_DIFF_INTRA:
			Correl_Type_String	=	"TWO FACTORS DIFF INTER DIFF INTRA";
			
			ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
			break;
	}

	if (TheX)
		delete	TheX;
	TheX	=	NULL;

	its_NbIter	=	nbiter;

	// -----------------------
	// STEP : GET RESULTS: and DISPLAY
	// -----------------------
	if (res)
		ICMLOG("BFGS - Compute CorrelationMatrix " << Correl_Type_String << " :: SUCCEEDED")
	else
		ICMLOG("BFGS - Compute_BetaBestFit_CorrelationMatrix " << Correl_Type_String << " :: FAILED")
	ICMLOG("BFGS::Results obtained after " << its_NbIter << " iterations and Val " << val)


}

// ---------------------------------------------------------------------------------------------------
// BETA BEST FIT
// ---------------------------------------------------------------------------------------------------

double	ICM_Correlation_Fitting::TheBetaBestFitFunction(double* X)
{
	int	i, j;
	double	Value;
	double	rho_i, rho_j;
	double	Sum;

	Sum	=	0.0;

	const ICM_QMatrix<double>& Matrix = its_CorrelationMatrix->GetMatrix();

//	vector<double>	TheParameters;
	ExtractRootsForBestFit(X);	//, TheParameters);

	for (i=1; i<its_CorrelationSize; i++)
	{
		rho_i	=	its_correl_parameters[i];

		for (j=0;j<i;j++)
		{
			rho_j	=	its_correl_parameters[j];
			Value	=	Matrix.Getvalue(i, j) - rho_i * rho_j;
			Value	*=	Value;
			Sum +=  Value;
		}
	}

	return	Sum;
}


void	ICM_Correlation_Fitting::TheBetaBestFitGradientFunction(double* X, double* Gradient)
{
	int	i;
	double	ResultUp, ResultDown, TheStep;
	
	for (i=0; i<its_dim; i++)
	{
		TheStep	=	its_Steps[i];
		// function evaluation
		X[i] += TheStep;		
		ResultUp	=	TheBetaBestFitFunction(X);

		// second function evluation
		X[i] -= 2.0 * TheStep;
		ResultDown	=	TheBetaBestFitFunction(X);

		// restore current point
		X[i] += TheStep;

		// set the gradient
		Gradient[i] = (ResultUp - ResultDown) / (2.0 * TheStep);
	}
}


// ---------------------------------------------------------------------------------------------------
// 2 FACTOR SAME INTER SAME INTRA BEST FIT
// ---------------------------------------------------------------------------------------------------

double	ICM_Correlation_Fitting::The2FSameInterSameIntraBestFitFunction(double* X)
{
	int	i, j;
	double	Value;
	double	Sum;

	int	iSecId;
	int	jSecId;

	Sum	=	0.0;

	const ICM_QMatrix<double>& Matrix = its_CorrelationMatrix->GetMatrix();

//	vector<double>	TheParameters;
	ExtractRootsForBestFit(X);	//, TheParameters);

	double	correl_inter;
	double	correl_intra;

	correl_inter	=	its_correl_parameters[0];
	correl_intra	=	its_correl_parameters[1];

	for (i=1; i<its_CorrelationSize; i++)
	{
		iSecId	=	GetSectorId(i);

		for (j=0;j<i;j++)
		{
			jSecId	=	GetSectorId(j);

			Value	=	Matrix.Getvalue(i, j);
			if (iSecId == jSecId)
				Value	-= correl_intra;
			else
				Value	-= correl_inter;

			Value	*=	Value;
			Sum +=  Value;
		}
	}

	return	Sum;
}


void	ICM_Correlation_Fitting::The2FSameInterSameIntraBestFitGradientFunction(double* X, double* Gradient)
{
	int	i;
	double	ResultUp, ResultDown, TheStep;
	
	for (i=0; i<its_dim; i++)
	{
		TheStep	=	its_Steps[i];
		// function evaluation
		X[i] += TheStep;		
		ResultUp	=	The2FSameInterSameIntraBestFitFunction(X);

		// second function evluation
		X[i] -= 2.0 * TheStep;
		ResultDown	=	The2FSameInterSameIntraBestFitFunction(X);

		// restore current point
		X[i] += TheStep;

		// set the gradient
		Gradient[i] = (ResultUp - ResultDown) / (2.0 * TheStep);
	}
}


void	ICM_Correlation_Fitting::ExtractRootsForBestFit(double* X) //, DoubleVector& Correl_Parameters)
{
	int	i;
//	its_correl_parameters.resize(its_dim);

	switch (its_Correlation_Fit_Type)
	{
		case	qCORREL_1F_SINGLE_VALUE:
			break;

		case	qCORREL_1F_BETAS:
			
			for (i=0; i<its_dim; i++)
				its_correl_parameters[i]	=	0.5 + atan(X[i]) / MY_PI;

			break;

		case	qCORREL_2F_SAME_INTER_SAME_INTRA:
			
			for (i=0; i<its_dim; i++)
				its_correl_parameters[i]	=	0.5 + atan(X[i]) / MY_PI;

			break;

		case	qCORREL_2F_SAME_INTER_DIFF_INTRA:

			ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
			break;

		case	qCORREL_2F_DIFF_INTER_DIFF_INTRA:
			
			ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
			break;
	}

}


void	ICM_Correlation_Fitting::Get_BestFit_Results(vector<double>& data)
{
	int	i;

	data.clear();
	data.resize(its_dim);

	switch (its_Correlation_Fit_Type)
	{
		case	qCORREL_1F_SINGLE_VALUE:
			data[0]	=	its_1F_Single_Value;
			break;

		case	qCORREL_1F_BETAS:
			
			for (i=0; i<its_dim; i++)
				data[i]	=	its_1F_Betas[i];

			break;

		case	qCORREL_2F_SAME_INTER_SAME_INTRA:

			data[0]	=	its_2F_Inter_Sector_Value;
			data[1]	=	its_2F_Intra_Sector_Value;

			break;

		case	qCORREL_2F_SAME_INTER_DIFF_INTRA:

			ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
			break;

		case	qCORREL_2F_DIFF_INTER_DIFF_INTRA:
			
			ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
			break;
	}

}

void	ICM_Correlation_Fitting::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else
	{
		fOut = ficOut;
	} 

	fprintf(fOut, "\n ======> Correlation Fitting:\n\n");
		

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}