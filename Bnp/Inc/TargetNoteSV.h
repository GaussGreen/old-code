#ifndef _TARGETNOTE_SV_H_
#define _TARGETNOTE_SV_H_

#include "TargetNote.h"

// ------------------------------------------------------------------------------------------------------------------- //
//
// Wrapper routine for pricing a TargetNote in LGMSV
//
char* TargetNoteMC_SV( 
						char *szUnd,
						TARN_Struct* tarn,
						TARN_AUX* aux,
						double **prod_val
				  );












#endif _TARGETNOTE_SV_H_
