/****************************************************************************/
/*        Memory allocation.	                                            */
/****************************************************************************/
/*        ALLOC.c                                                           */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "esl_alloc.h"
#include "esl_error.h"
#include "esl_macros.h"



/*****  MktVol_Init  **********************************************************/
/**
*       Initialize MKTVOL_DATA structure.
*/
void    MktVol_Init (
		MKTVOL_DATA    *mktvol_data /** Market Vol */
		)
{
	int    i;

	for (i=0; i<3; i++)
	{
        mktvol_data->Beta[i]  = 0e0;
        mktvol_data->Alpha[i] = 0e0;
        mktvol_data->Rho[i]   = 0e0;
    }

    mktvol_data->NbFactor          = 1;
    mktvol_data->FilterSpotVolFlag = FALSE;
    mktvol_data->VolUnit           = 1;
    mktvol_data->SmoothingFlag     = 'N';
    mktvol_data->TraceFlag         = 'Y';

    mktvol_data->ModelChoice       = FIX3_ORIGINAL;
    mktvol_data->IsNmrModel        = FALSE;
    mktvol_data->CalibSmileFlag    = TRUE;
    mktvol_data->NbNmr             = 0;
    mktvol_data->NbSigmaMQ         = -999;
    mktvol_data->NckMQ             = -999;
    
    return;

}  /* Dev_Init */



/*****  Opt_Out_Data_Init  **************************************************/
/**
 *		Initialize the OPT_OUT_DATA structure
 */
void     Opt_Out_Data_Init(
		OPT_OUT_DATA  *ood 
		)
{
	int  i,j;

	if (ood == NULL) return;

	ood->Option = 0.0;

	for (i=0; i<10; i++) ood->Price[i] = 0.0;

	for (i=0; i<ESL_NB_EVENT; i++) 
        {
            ood->prob_calc[i] = FALSE;

	    for (j=0; j<OPT_OUT_DATA_SIZE; j++)
            {
                 ood->prob_out_data[i].EventDate[j] = 0;
                 ood->prob_out_data[i].EventProb[j] = 0.;
            }
        }

	return;

} /* Opt_Out_Data_Init */

/*****  DR_Array  ***********************************************************/
/**
*       Allocation of memory for an array of arbitrary type. Returns void*.
*/
void    *DR_Array ( int     type    /** (I) Type         */ 
                   ,int     nl      /** (I) Lower bound  */
                   ,int     nh      /** (I) Higher bound */
		)
{
    switch (type)
    {
        case CRITDATE:
        {
            CRIT_DATE *v;

            v = (CRIT_DATE *) calloc ((unsigned) (nh-nl+1), sizeof (CRIT_DATE));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            /* Shift pointer and cast to void* */

            return ((void *) (v-nl));

        }         
        case INT:
        {
            int *v;

            v = (int *) calloc ((unsigned) (nh-nl+1), sizeof (int));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_PTR:
        {
            int **v;

            v = (int **) calloc ((unsigned) (nh-nl+1), sizeof (int *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_D_PTR:
        {
            int ***v;

            v = (int ***) calloc ((unsigned) (nh-nl+1), sizeof (int **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_T_PTR:
        {
            int ****v;

            v = (int ****) calloc ((unsigned) (nh-nl+1), sizeof (int ***));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG:
        {
            long *v;

            v = (long *) calloc ((unsigned) (nh-nl+1), sizeof (long));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG_PTR:
        {
            long **v;

            v = (long **) calloc ((unsigned) (nh-nl+1), sizeof (long *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG_D_PTR:
        {
            long ***v;

            v = (long ***) calloc ((unsigned) (nh-nl+1), sizeof (long **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE:
        {
            double *v;

            v = (double *) calloc ((unsigned) (nh-nl+1), sizeof (double));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE_PTR:
        {
            double **v;

            v = (double **) calloc ((unsigned) (nh-nl+1), sizeof (double *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE_D_PTR:
        {
            double ***v;

            v = (double ***) calloc ((unsigned) (nh-nl+1), sizeof (double **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }
        case CHAR:
        {
            char *v;

            v = (char *) calloc ((unsigned) (nh-nl+1), sizeof (char));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        } 
        case CHAR_PTR:
        {
            char **v;

            v = (char **) calloc ((unsigned) (nh-nl+1), sizeof (char *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        case TYPE_TPROB_0:
        {
            TPROB_0 *v;

            v = (TPROB_0 *) calloc ((unsigned) (nh-nl+1), sizeof (TPROB_0));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        case TYPE_TPROB_1:
        {
            TPROB_1 *v;

            v = (TPROB_1 *) calloc ((unsigned) (nh-nl+1), sizeof (TPROB_1));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        case TYPE_TPROB_2:
        {
            TPROB_2 *v;

            v = (TPROB_2 *) calloc ((unsigned) (nh-nl+1), sizeof (TPROB_2));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        case TYPE_SWAPRATE_DATA:
        {
            SWAPRATE_DATA *v;

            v = (SWAPRATE_DATA *) calloc ((unsigned) (nh-nl+1), sizeof (SWAPRATE_DATA));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        default:
        {
            DR_Error ("Program bug: unsupported type in DR_Array!");
            return (NULL);
        }                          
    }  /* switch */

}  /* DR_Array */

/*****  DR_Matrix  **********************************************************/
/**
 *      Allocation of memory for a matrix of arbitrary type. Returns void *.
 */
void  *DR_Matrix (int     type    /** (I) Type                */
                 ,int     nrl     /** (I) Lower row index     */
                 ,int     nrh     /** (I) Higher row index    */
                 ,int     ncl     /** (I) Lower column index  */
                 ,int     nch     /** (I) Higher column index */
		)
{
    void  **Matx = NULL;
    int     i;
    int     type_ptr;

    switch (type)
    {
    case INT:
        type_ptr = INT_PTR;
        break;
    case LONG:
        type_ptr = LONG_PTR;
        break;
    case DOUBLE:
        type_ptr = DOUBLE_PTR;
        break;
    default:
        DR_Error ("DR_Matrix: unsupported type.");
        goto RETURN;
    }  /* switch */

    Matx = (void **) DR_Array (type_ptr, nrl, nrh);

    if (Matx == NULL)
    {
        DR_Error ("DR_Matrix: allocation failure #1!");
        goto RETURN;
    }

    for (i = nrl; i <= nrh; i++)
    {
        Matx[i] = DR_Array (type, ncl, nch);

        if (Matx[i] == NULL)
        {
            DR_Error ("DR_Matrix: allocation failure #2!");
            goto RETURN;
        }
    }

    return ((void *) Matx);

RETURN:

    /* allocation failed, release memory alloc'd so far */
    Free_DR_Matrix (Matx, type, nrl, nrh, ncl, nch);

    return (NULL);

} /* DR_Matrix */


/*****  DR_Cube  ****************************************************************/
/**
 *       Allocation of memory for a 3D array of arbitrary type. Returns void *.
 */
void    *DR_Cube (  int     type    /** (I) Type                */
                   ,int     nl      /** (I) Lower bound         */
                   ,int     nh      /** (I) Higher bound        */
                   ,int     nrl     /** (I) Lower row index     */
                   ,int     nrh     /** (I) Higher row index    */
                   ,int     ncl     /** (I) Lower column index  */
                   ,int     nch     /** (I) Higher column index */
		)
{
   void ***Cube = NULL;
    int     i, j;
    int     type_ptr, type_dbl_ptr;

    switch (type)
    {
    case INT:
        type_dbl_ptr = INT_D_PTR;
        type_ptr     = INT_PTR;
        break;
    case LONG:
        type_dbl_ptr = LONG_D_PTR;
        type_ptr     = LONG_PTR;
        break;
    case DOUBLE:
        type_dbl_ptr = DOUBLE_D_PTR;
        type_ptr     = DOUBLE_PTR;
        break;
    default:
        DR_Error ("DR_Cube: unsupported type.");
        goto RETURN;
    }  /* switch */

    Cube = (void ***) DR_Array (type_dbl_ptr, nl, nh);

    if (Cube == NULL)
    {
        DR_Error ("DR_Cube: allocation failure #1!");
        goto RETURN;
    }

    for (i = nl; i <= nh; i++)
    {
        Cube[i] = (void **) DR_Array (type_ptr, nrl, nrh);
        if (Cube[i] == NULL)
        {
            DR_Error ("DR_Cube: allocation failure #2!");
            goto RETURN;
        }
        
        for (j = nrl; j <= nrh; j++)
        {
            Cube[i][j] = DR_Array (type, ncl, nch);
            
            if (Cube[i][j] == NULL)
            {
                DR_Error ("DR_Cube: allocation failure #3!");
                goto RETURN;
            }
        } /* j */
    } /* i */

    return ((void *) Cube);

RETURN:

    /* allocation failed, release memory alloc'd so far */
    Free_DR_Cube (Cube, type, nl, nh, nrl, nrh, ncl, nch);

    return (NULL);

}  /* DR_Cube */


/*****  Free_DR_Array  ******************************************************/
/**
*       Free DR array.
*/
int     Free_DR_Array ( void   *Array   /** (I) Array        */
                       ,int     type    /** (I) Type         */
                       ,int     nl      /** (I) Lower bound  */
                       ,int     nh      /** (I) Higher bound */
		)
{
    /* To avoid warning message */
    nh += 0;

    switch (type)
    {
        case CRITDATE:
        {
            CRIT_DATE *v = (CRIT_DATE *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }  
        case INT:
        {
            int *v = (int *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_PTR:
        {
            int **v = (int **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_D_PTR:
        {
            int ***v = (int ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_T_PTR:
        {
            int ****v = (int ****) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG:
        {
            long *v = (long *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG_PTR:
        {
            long **v = (long **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG_D_PTR:
        {
            long ***v = (long ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE:
        {
            double *v = (double *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE_PTR:
        {
            double **v = (double **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE_D_PTR:
        {
            double ***v = (double ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }
        case CHAR:
        {
            char *v = (char *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case CHAR_PTR:
        {
            char **v = (char **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }  
        case TYPE_TPROB_0:
        {
            TPROB_0 *v = (TPROB_0 *) Array;

            if (v != NULL)
            {
                free ((TPROB_0 *) (v+nl));
            }

            break;
        }                          
        case TYPE_TPROB_1:
        {
            TPROB_1 *v = (TPROB_1 *) Array;

            if (v != NULL)
            {
                free ((TPROB_1 *) (v+nl));
            }

            break;
        }                          
        case TYPE_TPROB_2:
        {
            TPROB_2 *v = (TPROB_2 *) Array;

            if (v != NULL)
            {
                free ((TPROB_2 *) (v+nl));
            }

            break;
        }                          
        case TYPE_SWAPRATE_DATA:
        {
            SWAPRATE_DATA *v = (SWAPRATE_DATA *) Array;

            if (v != NULL)
            {
                free ((SWAPRATE_DATA *) (v+nl));
            }

            break;
        }  
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Array!");
            return (FAILURE);
        }                          
    }  /* switch */

    Array = NULL;

    return (SUCCESS);

}  /* Free_DR_Array */



/*****  Free_DR_Matrix  *****************************************************/
/**
*       Free DR matrix.
*/
int     Free_DR_Matrix (void   *Matrix      /** (I) Matrix              */
                       ,int     type        /** (I) Type                */
                       ,int     nrl         /** (I) Lower row index     */
                       ,int     nrh         /** (I) Higher row index    */
                       ,int     ncl         /** (I) Lower column index  */
                       ,int     nch         /** (I) Higher column index */
		)
{
    int i;

    /* To avoid warning message */
    nch += 0;

    switch (type)
    {
        case INT:
        {
            int **m = (int **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        case LONG:
        {
            long **m = (long **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        case DOUBLE:
        {
            double **m = (double **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Matrix!");
            return (FAILURE);
        }                          
    }  /* switch */

    Matrix = NULL;

    return (SUCCESS);

}  /* Free_DR_Matrix */



/*****  Free_DR_Cube  *******************************************************/
/**
*       Free 3D DR array.
*/
int     Free_DR_Cube (  void   *Cube    /** (I) Cube                */
                       ,int     type    /** (I) Type                */
                       ,int     nl      /** (I) Lower bound         */
                       ,int     nh      /** (I) Higher bound        */
                       ,int     nrl     /** (I) Lower row index     */
                       ,int     nrh     /** (I) Higher row index    */
                       ,int     ncl     /** (I) Lower column index  */
                       ,int     nch     /** (I) Higher column index */
		)
{
    int i;
    int j;

    /* To avoid warning message */
    nch += 0;

    switch (type)
    {
        case INT:
        {
            int ***c = (int ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                          
        case LONG:
        {
            long ***c = (long ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                          
        case DOUBLE:
        {
            double ***c = (double ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                        
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Cube!");
            return (FAILURE);
        }                          
    }  /* switch */

    Cube = NULL;

    return (SUCCESS);

}  /* Free_DR_Cube */

void InitZeroCurve(T_CURVE* crv)
{
#ifdef ESL_NEW_CURVE
    irxZeroCurveInit(crv);
#endif
}

void DestroyZeroCurve(T_CURVE* crv)
{
#ifdef ESL_NEW_CURVE
    irxZeroCurveDestroy(crv);
#endif
}


