/***
	calibrate.h

	SORT 
	
	Author Adam H. Litke
        Uses utilities by E. Auld, G. Amblard, and A. Litke

	Some routines written by K L Chau

	Two Factor Routines written by OVE
***/
Err fit(SrtUndPtr und,
        SrtMdlPtr model,
        SwapDP *sdp,
        double *strk,
        double *bnd_strk,
        int 	*t,
        SrtReceiverType	*rp,
        double *price,
        double *vega,
        int    num_inp,
        double **v,
        double **l,
        int    *len,
        FILE   *f,
        double *answer,
        long *volDateOut,
	double *volOut,
	long *tauDateOut,
	double *tauOut,
	int frz_tau); 
                                      
Err tf_global_fit(
		SrtUndPtr und,
        SrtMdlPtr model,
        SwapDP *sdp,
        double *strk,
        double *bnd_strk,
        int    *t,
		SrtReceiverType *rec_pay,
        double *price,
        double *vega,
        int    num_inp,
        double **v,
        double **l,
        int    *len,
        FILE   *f,
        double *answer) ; 

Err tf_const_cor_fit(
		SrtUndPtr und,
        SrtMdlPtr model,
        SwapDP *sdp,
        double *strk,
        double *bnd_strk,
        int    *t,
		SrtReceiverType *rec_pay,
        double *price,
        double *vega,
        int    num_inp,
        double **v,
        double **l,
        int    *len,
        FILE   *f,
        double *answer) ; 

Err tf_fix_cor_fit(
		SrtUndPtr und,
        SrtMdlPtr model,
        SwapDP *sdp,
        double *strk,
        double *bnd_strk,
        int    *t,
		SrtReceiverType *rec_pay,
        double *price,
        double *vega,
        int    num_inp,
        double **v,
        double **l,
        int    *len,
        FILE   *f,
        double *answer) ; 

Err fit_fn(
		SrtUndPtr          und,
		SrtMdlPtr          model,
		SwapDP             *sdp,
		double             *strk ,
		double             *bnd_strk,
		int                *t,
		SrtReceiverType    *rp,
		double             *price,
		double             *vega,
		String             *refratecode,
		int                num_inp,
		double             **vol,
		int                nvols,
		double             **l,
		int                ntaus,
		FILE               *f,
		double             *answer,
		void               *func,
		int                nfunc,
		SRT_Boolean            flag,
		int 	           freeze_tau);


Err tf_global_fit_basis_fn(
		SrtUndPtr und,
        SrtMdlPtr model,
        SwapDP	*sdp,
        double	*strk,
        double	*bnd_strk,
        int	*t,
		SrtReceiverType *rec_pay,
        double	*price,
        double	*vega,
		String  *refratecode,
        int 	num_inp,
        double	**vol,
        double	**l,
        int	*len,
        FILE	*f,
        double	*answer,
        int     **funcs, 
	int 	nfunc[3],
	int 	freeze_tau,
	SRT_Boolean use_predefined_fn,
	double  **user_fn_range,
	int	user_fn_length);




double basis_fn_coeff(int index);
double tf_basis_fn_coeff(int fitted_set, int index);

#define NUM_OF_BASIS_FUNC 	7
#define NUM_FIT_PARAM     	4 

                                      
