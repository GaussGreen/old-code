#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123.h"
#include "fix123head.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Fix3_MC_ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

void Fix3_MC_StartTreeSim(TREESIM_DATA *ts)
{
    int p;

    /* all paths start at node (0, 0, 0) */

    for (p=0; p<ts->NbPathAll; ++p)
    {
        ts->i[p] = 0;
        ts->j[p] = 0;
        ts->k[p] = 0;
    }

    ts->xx = -1 - ts->seed; /* input seed is nonnegative;     */
                            /* we need a negative actual seed */
}

int Fix3_MC_UpdateTreeSim(int                dim,
                          int                t,
                          TREESIM_DATA      *ts,
                          FIX3_DEV_DATA     *dev_data,
                          FIX3_TREE_DATA    *tree_data)
{
    int iold, jold, kold;
    int inew, jnew, knew;
    int z;

    int offset1D;
    int offset2D;
    int offset3D;

    int *Shift1;
    int *Shift2;
    int *Shift3;
    
    int     Top1,   Bottom1;
    int    *Top2,  *Bottom2;
    int   **Top3, **Bottom3;

    /* not used
    TPROB *p;
    TPROB *q;
    TPROB *r;
    */

    double *p0L;
    double *pdL;

    double p0;
    double pd;

    double q0;
    double qd;

    double r0;
    double rd;

    double y;

    int status = FAILURE;

    if (t >= tree_data->NbTP)
    {
        return SUCCESS;
    }

    if ((dim != 1) &&
        (dim != 2) &&
        (dim != 3))
    {
        DR_Error("Tree dim. must be 1,2, or 3. (UpdatePathIndexes)\n");
        goto done;
    }

    offset1D = Fix3_Node_Offset(1, 0, 0, t, tree_data);

    Shift1  = dev_data->Shift1 + offset1D;

    if (dim == 1)
    {
        pdL = dev_data->pd + offset1D;
        p0L = dev_data->p0 + offset1D;
    }

    Top1 = tree_data->Top1[t+1];
    Bottom1 = tree_data->Bottom1[t+1];

    Top2 = tree_data->Top2[t+1];
    Bottom2 = tree_data->Bottom2[t+1];

    Top3 = tree_data->Top3[t+1];
    Bottom3 = tree_data->Bottom3[t+1];

    for (z=0; z<ts->NbPathAll; ++z)
    {
        iold = ts->i[z];

        inew = iold + Shift1[iold];

        if (dim == 1) 
        {
            pd = pdL[iold];
            p0 = p0L[iold];
        }
        else
        {
            int offset;
            double q0u, q00, q0d, qdu, qd0, qdd;

            offset2D = Fix3_Node_Offset(2, iold, 0, t, tree_data);

            Shift2 = dev_data->Shift2 + offset2D;

            jold = ts->j[z];

            jnew = jold + Shift2[jold];

            offset = offset2D +jold;

            q0u = *(dev_data->q0u + offset);
            q00 = *(dev_data->q00 + offset);
            q0d = *(dev_data->q0d + offset);
            qdu = *(dev_data->qdu + offset);
            qd0 = *(dev_data->qd0 + offset);
            qdd = *(dev_data->qdd + offset);

            /* Back out pd, p0, qd, q0 */
            pd =  qdu + qd0 + qdd;
            p0 =  q0u + q00 + q0d;

            qd = qdd / pd;
            q0 = qd0 / pd;
        }


        y = Fix3_MC_ran2(&(ts->xx));

        if (y < pd)
        {
            inew -= 1;
        }
        else if (y > pd + p0)
        {
            inew += 1;
        }

        inew = MINMAX(inew, Bottom1, Top1);

        ts->i[z] = inew;

        if (dim > 1)
        {
            y = Fix3_MC_ran2(&(ts->xx)) * pd;

            if (y < qd)
            {
                jnew -= 1;
            }
            else if (y > qd + q0)
            {
                jnew += 1;
            }

            jnew = MINMAX(jnew, Bottom2[inew], Top2[inew]);

            ts->j[z] = jnew;

            if (dim > 2)
            {
                offset3D = Fix3_Node_Offset(3, iold, jold, t, tree_data);

                Shift3 = dev_data->Shift3 + offset3D;

                kold = ts->k[z];

                knew = kold + Shift3[kold];

                rd = *(dev_data->rd + offset3D + kold);
                r0 = *(dev_data->r0 + offset3D + kold);

                y = Fix3_MC_ran2(&(ts->xx));

                if (y < rd)
                {
                    knew -= 1;
                }
                else if (y > rd + r0)
                {
                    knew += 1;
                }

                knew = MINMAX(knew, Bottom3[inew][jnew], Top3[inew][jnew]);

                ts->k[z] = knew;
            }
        }
    }

    status = SUCCESS;

done:

    return status;
}

int IsGreater(const void *a, const void *b)
{
    double x = *(const double *)a - *(const double *)b;
    if (x<0) return -1;
    else if (x>0) return 1;
    else return 0;
}

int Fix3_MC_ProcessSamples(TREESIM_DATA *ts)
{
    int i;
    int j;
    int k;

    int NbPathAll = ts->NbPathAll;
    int NbPathSub = ts->NbPathSub;

    int NbStateLeft;

    /* our sample will be:                              */
    /* 1) points evenly spaced in (sorted) index space, */
    /*    including endpoints                           */
    /* 2) points evenly spaced in value space           */

    int NbPathSubHalf = NbPathSub / 2;

    double index;

    double a;
    double b;

    double x = (double)(NbPathAll - 1) / (double)(NbPathSubHalf - 1);
    double y;

    double value;

    for(i=0; i<ts->NbPathDate; ++i)
    {
        /* sort full sample */

        qsort((void *)ts->PathAll[i], ts->NbPathAll, sizeof(double), &IsGreater);

        /* select subsample, including endpoints */

        index = 0;
        j     = 0;
        k     = 0;

        while (j < ts->NbPathAll)
        {
            ts->PathSub[i][k] = ts->PathAll[i][j];

            k     += 1;
            index += x;

            j = NEAR_INT(index);
        }

        /* eliminate duplicates */

        j = 0;
        k = 0;

        ts->State[i][j] = ts->PathSub[i][k];
    
        while (k < NbPathSubHalf)
        {
            while ((k < NbPathSubHalf) &&
                   (IS_EQUAL(ts->State[i][j], ts->PathSub[i][k])))
            {
                ++k;
            }

            if (k < NbPathSubHalf)
            {
                ++j;
                ts->State[i][j] = ts->PathSub[i][k];
            }
        }

        ts->NbState[i] = j+1;

        /* add points evenly spaced between StateMin and StateMax */

        NbStateLeft = ts->NbPathSub - ts->NbState[i];

        a = ts->StateMin[i];
        b = ts->StateMax[i];

        y = (b - a) / (double)(NbStateLeft - 1);

        value = a;

        for (k=ts->NbState[i]; k<NbPathSub; ++k)
        {
            ts->State[i][k] = value;

            value += y;
        }

        /* sort the final sample */

        qsort((void *)ts->State[i], NbPathSub, sizeof(double), &IsGreater);

        /* eliminate duplicates again, this time "in situ" */

        j = 0;
        k = 0;

        ts->State[i][j] = ts->State[i][k];
    
        while (k < NbPathSub)
        {
            while ((k < NbPathSub) &&
                   (IS_EQUAL(ts->State[i][j], ts->State[i][k])))
            {
                ++k;
            }

            if (k < NbPathSub)
            {
                ++j;
                ts->State[i][j] = ts->State[i][k];
            }
        }

        ts->NbState[i] = j+1;

        ts->MaxNbState = MAX(ts->MaxNbState, ts->NbState[i]);
    }

    return SUCCESS;
}

int Fix3_MC_CheckSamples(TREESIM_DATA *ts)
{
    int i;

    for (i=0; i<ts->NbPathDate; ++i)
    {
        if ((ts->StateMin[i] > ts->PathAll[i][0] + TINY) ||
            (ts->StateMax[i] < ts->PathAll[i][ts->NbPathAll-1] - TINY))
        {
            DR_Error("Wide state variable limits at %ld\n"
                     "fail to encompass MC sample.\n", ts->PathDate[i]);

            return FAILURE;
        }
    }

    return SUCCESS;
}

int Fix3_MC_CheckDate(FIX3_TREE_DATA *tree_data,
                      int             t,
                      TREESIM_DATA   *ts,
                      int             SrIdx,
                      int             SrDateIdx)
{
    int status = FAILURE;

    long CurrTreeDate, CurrMcDate;

    if ((SrIdx <  0) ||
        (SrIdx >= ts->NbSwapRateData)) goto done;

    if ((SrDateIdx <  0) ||
        (SrDateIdx >= ts->SwapRateData[SrIdx].NbResetDate)) goto done;

    CurrTreeDate = tree_data->TPDate[t];
    CurrMcDate   = ts->SwapRateData[SrIdx].ResetDate[SrDateIdx];

    if (CurrTreeDate != CurrMcDate) goto done;

    status = SUCCESS;

done:

    return status;
}

int Fix3_MC_GetSliceValue(double         *Val,
                          FIX3_TREE_DATA *tree_data,
                          TREESIM_DATA   *ts,
                          TSLICE          Slice,
                          int             dim,
                          int             TrDateIdx,
                          int             path)
{
    int status = FAILURE;
    int arg;

    TSLICE SliceL;

    SliceL = Slice +

             Fix3_Node_Offset(dim,
                              ts->i[path],
                              ts->j[path],
                              TrDateIdx,
                              tree_data);
 
    switch(dim)
    {
        case 1: arg = ts->i[path]; break;
        case 2: arg = ts->j[path]; break;
        case 3: arg = ts->k[path]; break;
        default:
        {
            DR_Error("Coding error: could not retrieve slice value.\n");
            goto done;
        }
    }

    *Val = SliceL[arg];

    status = SUCCESS;

done:

    return status;
}
