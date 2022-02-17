#ifndef __SRTGRFNMAIN3FQUANTO_H
#define __SRTGRFNMAIN3FQUANTO_H

char *SrtGrfn3FQuantoTree(
						  char*		und3fquanto,
						  double*	correl_dates,
						  int		n_correl_dates,
						  double***	correl_matrix_5x5_ts,
						  int		numeventdates, 
						  long*		eventdates,
						  long		tableauRows,
						  long		tableauCols,
						  char***	tableauStrings,
						  int**		tableauMask,
						  long		auxWidth,
						  long*		auxLen,
						  double**	aux,
						  long		num_stp, 
						  int*		num_prod, 
						  int		discount,
						  double**	prod_val);

#endif
