#ifndef ESL_ALLOC_DOT_H
#define ESL_ALLOC_DOT_H

/** NOTE: This file should be only included through 'esl_alloc.c'
 */


#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"

#ifdef  __cplusplus
extern "C" {
#endif


/*****  DR_Array  ***********************************************************/
/**
*       Allocation of memory for an array of arbitrary type. Returns void*.
*/
void    *DR_Array ( int     type    /** (I) Type         */ 
                   ,int     nl      /** (I) Lower bound  */
                   ,int     nh      /** (I) Higher bound */
		);

/*****  DR_Matrix  **********************************************************/
/**
 *      Allocation of memory for a matrix of arbitrary type. Returns void *.
 */
void  *DR_Matrix (int     type    /** (I) Type                */
                 ,int     nrl     /** (I) Lower row index     */
                 ,int     nrh     /** (I) Higher row index    */
                 ,int     ncl     /** (I) Lower column index  */
                 ,int     nch     /** (I) Higher column index */
		);

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
		);

/*****  Free_DR_Array  ******************************************************/
/**
*       Free DR array.
*/
int     Free_DR_Array ( void   *Array_   /** (I) Array        */
                       ,int     type    /** (I) Type         */
                       ,int     nl      /** (I) Lower bound  */
                       ,int     nh      /** (I) Higher bound */
		);

/*****  Free_DR_Matrix  *****************************************************/
/**
*       Free DR matrix.
*/
int     Free_DR_Matrix (void   *Matrix_      /** (I) Matrix              */
                       ,int     type        /** (I) Type                */
                       ,int     nrl         /** (I) Lower row index     */
                       ,int     nrh         /** (I) Higher row index    */
                       ,int     ncl         /** (I) Lower column index  */
                       ,int     nch         /** (I) Higher column index */
		);

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
		);


void InitZeroCurve(T_CURVE* crv);
void DestroyZeroCurve(T_CURVE* crv);


#ifdef  __cplusplus
}
#endif


#endif


