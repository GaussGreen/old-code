/* -------------------------------------------------------------------------
   FILENAME		: num_h_vector_tools.h
   AUTHOR		: Toto 16Dec1999

   ------------------------------------------------------------------------- */

#ifndef NUM_H_VECTOR_TOOLS_H
#define NUM_H_VECTOR_TOOLS_H

/*	Fill evenly a vector so that it contains a specified number of points */
void num_f_fill_vector(
    /*	Size of the input vector	*/
    int p,
    /*	Vector byref	*/
    double** v,
    /*	Required size of the output vector	*/
    int n,
    /*	Required size of the first section	*/
    int n0);

/*	Another algorithm... */
void num_f_fill_vector_newalgo(
    /*	Initial size of vector	*/
    int* p,
    /*	Vector byref	*/
    double** v,
    /*	Required size of the output vector	*/
    int n);

/*	Fill evenly a vector so that the space between two consecutive points is less
        than a specified value */
void num_f_fill_vector_maxtime(
    /*	Initial size of vector	*/
    int* p,
    /*	Vector byref	*/
    double** v,
    /*	Required maxtime	*/
    double maxtime);

/*	Add a double	*/
void num_f_add_number(
    /*	Initial size of vector	*/
    int* p,
    /*	Initial vector byref	*/
    double** v1,
    /*	Number to be addded	*/
    double d);

/*	Concatenate two vectors	*/
void num_f_concat_vector(
    /*	Initial size of vector	*/
    int* p,
    /*	Initial vector byref	*/
    double** v1,
    /*	Size of vector to be concatenated	*/
    int n,
    /*	Vector to be concatenated	*/
    double* v2);

/*	Sort vector	*/
void num_f_sort_vector(
    /*	Size of vector	*/
    int p,
    /*	Vector	*/
    double* v);

/*	Remove identical values	*/
void num_f_unique_vector(
    /*	Size of vector byref	*/
    int* p,
    /*	Vector	*/
    double* v);

void num_f_sort_matrix(
    /*	Size of vector			*/
    int p,
    /*	Number of column		*/
    int n,
    /*	Column on which we sort	*/
    int col,
    /*	Vector	*/
    double** v);

#endif
