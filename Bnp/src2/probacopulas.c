/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "math.h"
#include "num_h_proba.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: num_f_Gaussian_copulas(...)
 *
 * DESCRIPTION  	: Returns the correlation matrix used for
 *                 gaussian copula thanks to the Sperman depedency factor
 *
 * CALLS			: None
 *
 *
 * PARAMETERS
 *	INPUT		: spot			- spot underlying price
 *				: barrier      	- barrier level to stay under
 *				: vol         	- annual volatility
 *				: mat         	- maturity, in years
 *				: disc        	- discount factor for the domestic currency
 *
 * RETURNS		: proba			- probability
 *
 *******************************************************************************/

static int CompFunc(const void* vp, const void* vq)
{
    const double* p = vp;
    const double* q = vq;

    if ((*p) > (*q))
        return (+1);
    else if ((*p) < (*q))
        return (-1);
    else
        return (0);
}

/**********************************************************************************

   In order to estimate the empirical Copula we compute

                                     number of pairs (x,y) in the sample such as x<=x(i) and y<=y(i)
                                Cn = ---------------------------------------------------------------
                                                                                                n

***********************************************************************************/

double EmpiricalCopula(
    double** InitialDatas,
    double   SortedData1,
    double   SortedData2,
    long     index1,
    long     index2,
    long     nDates)
{
    long   i;
    double result;

    result = 0.0;
    for (i = 0; i < nDates; i++)
    {
        if ((InitialDatas[index1][i] <= SortedData1) && (InitialDatas[index2][i] <= SortedData2))
            result += 1.0;
    }

    return result / nDates;
}

/**********************************************************************************

   In order to estimate the empirical Copula frequency we compute

                                     1/n if (x(i),y(i)) is an element of the sample
                                cn =
                                         0 otherwise

  Remark : One could notice that Cn(i/n,i/j) = sum{p,q} cn(p/n,q/n)

***********************************************************************************/

double EmpiricalCopulaFrequency(
    double** InitialDatas,
    double   SortedData1,
    double   SortedData2,
    long     index1,
    long     index2,
    long     nDates)
{
    long   i;
    double result;

    result = 0.0;
    i      = 0;
    while ((i < nDates) && (result == 0))
    {
        if ((InitialDatas[index1][i] == SortedData1) && (InitialDatas[index2][i] == SortedData2))
            result = 1;

        i++;
    }

    return result / nDates;
}

double SpearmansRho(
    double** InitialDatas, double** SortedDatas, long index1, long index2, long nDates)
{
    long   i, j;
    double result;
    double EmpCopul;

    result = 0.0;

    for (i = 0; i < nDates; i++)
    {
        for (j = 0; j < nDates; j++)
        {
            EmpCopul = EmpiricalCopula(
                InitialDatas,
                SortedDatas[index1][i],
                SortedDatas[index2][j],
                index1,
                index2,
                nDates);

            result += EmpCopul - (i + 1.0) * (j + 1.0) / nDates / nDates;
        }
    }

    return result * 12.0 / (nDates * nDates - 1.0);
}

Err num_f_Gaussian_copulas(
    double**  Datas, /* Datas[0..nDim - 1][0..nNumDates - 1] */
    long      nDim,
    long      nNumDates,
    double*** CopulaCorrelation) /* no allocation inside !
                                   assumed that it is [0..nDim][0..nDim] */
{
    Err      err = NULL;
    long     i, j;
    double** SortedDataAndIndex  = NULL;
    double** InitialDataAndIndex = NULL;
    double   dTempCovar;

    /* Inintialisation of the Array of AtomData's */
    InitialDataAndIndex = (double**)malloc(nDim * sizeof(double*));
    SortedDataAndIndex  = (double**)malloc(nDim * sizeof(double*));

    for (i = 0; i < nDim; i++)
    {
        InitialDataAndIndex[i] = (double*)malloc(nNumDates * sizeof(double));
        SortedDataAndIndex[i]  = (double*)malloc(nNumDates * sizeof(double));
        for (j = 0; j < nNumDates; j++)
        {
            InitialDataAndIndex[i][j] = Datas[i][j];
            SortedDataAndIndex[i][j]  = Datas[i][j];
        }
    }

    /* Sort the Datas */
    for (i = 0; i < nDim; i++)
        qsort(SortedDataAndIndex[i], nNumDates, sizeof(double), CompFunc);

    for (i = 0; i < nDim; i++)
    {
        (*CopulaCorrelation)[i][i] = 1;
        for (j = i + 1; j < nDim; j++)
        {
            dTempCovar = SpearmansRho(InitialDataAndIndex, SortedDataAndIndex, i, j, nNumDates);

            /* Get the Correlation for the Gaussian Copulas */
            (*CopulaCorrelation)[i][j] = 2.0 * sin(SRT_PI * dTempCovar / 6.0);
            (*CopulaCorrelation)[j][i] = (*CopulaCorrelation)[i][j];
        }
    }

    /* Free the memory */
    for (i = 0; i < nDim; i++)
        srt_free(InitialDataAndIndex[i]);
    srt_free(InitialDataAndIndex);

    for (i = 0; i < nDim; i++)
        srt_free(SortedDataAndIndex[i]);
    srt_free(SortedDataAndIndex);

    return err;
}
