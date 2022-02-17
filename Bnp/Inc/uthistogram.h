/* -------------------------------------------------------------------------------------
   
   FILE NAME : 	  UTHISTOGRAM.c	 			

   PURPOSE:       functions to store data for and create histograms

   ------------------------------------------------------------------------------------- */

#ifndef UTHISTOGRAM_H
#define UTHISTOGRAM_H

/* The linked list to store the histograms*/
typedef SrtListHdr SrtHistList, *SrtHistListPtr;
/*An element of the linked list of histogram*/
typedef SrtListAtom SrtHistAtom; 

/* The Data structure*/
typedef struct
{
	long num_path;
	double *data;
}
SrtHistData;

/* The built HISTOGRAM structure */
typedef struct
{
	char name[32];
	int *dens;
	double *abs;
	long num_seg;
}
SrtHistStr, *SrtHistStrPtr ;

/* --------------------------------  The Histogram Type -------------------------------- */

typedef enum srthistotype{
	SCATTERED,
	HISTOGRAM
} SrtHistoType;

/* ------------------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------------
                   SIX FUNCTIONS TO BE USED FOR AN EXTERNAL PROGRAM 
   ------------------------------------------------------------------------------------- */

/* ---- For the Histograms management ---- */

/* Function to be called when loading the executable */
Err destroy_histogram_list();

Err create_histogram_list( String histo_list_name);
	

Err srt_f_ininewhist(long num_path);

Err srt_f_sendhistval(String name, long path, double value);

int srt_f_hislsthas(String name);


/* - Creates the effective histogram (X or X-Y)  fromt he data stored - */ 
Err srt_f_makehistogram(
		String name_x,
		long   num_seg_x,
		String name_y,
		long   num_seg_y,
		/* OUTPUT */
		double **histo_x,
		long   *x_size,
		double **histo_y,
		long   *y_size,
		double ***histo_values);

/*
In histo_abs(0), we store the AVERAGE of data 
	and in histo_abs(-1) the STD_DEV of data
In histo_dens(0), we store the average_density for the graph=max(density)
	and in histo_dens(-1) the std_dev_density = 0.5 * max(density)
*/


/* -------------------------- Internal Functions --------------------------------------- */

/** Puts the input SrtHistData  into the _srt_global_histlistptr list	
and reference it as name**/
Err srt_f_addhistdata(String name,SrtHistData *data);

/** Picks up in _srt_global_histlistptr the data_base necessary 
to create the histogram referenced with name**/ 
SrtHistData *srt_f_gethistdata(String name);

/** Allocates space for the referenced SrtHistData and the data inside
and adds this SrtHistData as a SrtObject into _srt_global_histlistptr
under this new name  **/
SrtHistData *srt_f_crehistdata(String name, long npath);

/** Free the space used by SrtHistData and the data inside**/
Err srt_f_frehistdata(SrtHistData *data);

/** Same but for use as a freee function in a SrtList */
Err srt_f_frehistval(void *histdata);

/** Allocates space for a *SrtHistStr **/
Err srt_f_crehiststr(SrtHistStr **histo);

/** Free  contents and addresses of a *SrtHistStr **/
Err srt_f_frehiststr(SrtHistStr *histo);




/** Allocates space for  a *SrtHistList **/
Err srt_f_crehistlist(SrtHistList **hist_list, String hist_name);

/** Free contents (SrtHistData) and addresses of a SrtHistList	**/
Err srt_f_frehistlist(SrtHistListPtr *hist_list);


#endif

