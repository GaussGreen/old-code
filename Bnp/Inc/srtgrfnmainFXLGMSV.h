#ifndef SRTGRFNMAINFXLGMSV_H
#define	SRTGRFNMAINFXLGMSV_H


#include "srt_h_all.h"
#include "grf_h_all.h"

char *SrtGrfnFXLGMSVMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					// for Optimisation of exercise boundary 
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnQTOLGMSVMC(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnQTOLGMSV2FMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					// for Optimisation of exercise boundary 
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnQTOLGMSV1FMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					// for Optimisation of exercise boundary 
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

#endif
