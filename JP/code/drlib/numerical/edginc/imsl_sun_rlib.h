/*  Use single precision functions on sun, where possible */



static float arg1_f_, arg2_f_;

#define ARG1(X_F_)		(arg1_f_=(X_F_),&arg1_f_)

#define ARG2(X_F_)		(arg2_f_=(X_F_),&arg2_f_)



#ifdef sparc

static union {double d_f_;  float r_f_} u_f_;

#define ANS(Y)		(u_f_.d_f_=(Y),u_f_.r_f_)

#endif



#ifdef mc68020

static float val_f_;

#define ANS(Y)	    (*(int *)(&val_f_)=(Y),val_f_)

#endif



#define FFUN(F,X1_F_)		ANS(F(ARG1(X1_F_)))

#define FFUN2(F,X1_F_,X2_F_)	ANS(F(ARG1(X1_F_),ARG2(X2_F_)))



#define sin(X_F_)	FFUN(r_sin_,X_F_)

#define fabs(X_F_)	FFUN(r_fabs_,X_F_)

#define rint(X_F_)	FFUN(r_rint_,X_F_)

#define hypot(X_F_)	FFUN(r_hypot_,X_F_)

#define sqrt(X_F_)	FFUN(r_sqrt_,X_F_)

#define asinh(X_F_)	FFUN(r_asinh_,X_F_)

#define acosh(X_F_)	FFUN(r_acosh_,X_F_)

#define atanh(X_F_)	FFUN(r_atanh_,X_F_)

#define exp(X_F_)	FFUN(r_exp_,X_F_)

#define log(X_F_)	FFUN(r_log_,X_F_)

#define log10(X_F_)	FFUN(r_log10_,X_F_)

#define sin(X_F_)	FFUN(r_sin_,X_F_)

#define cos(X_F_)	FFUN(r_cos_,X_F_)

#define tan(X_F_)	FFUN(r_tan_,X_F_)

#define asin(X_F_)	FFUN(r_asin_,X_F_)

#define acos(X_F_)	FFUN(r_acos_,X_F_)

#define atan(X_F_)	FFUN(r_atan_,X_F_)

#define sinh(X_F_)	FFUN(r_sinh_,X_F_)

#define cosh(X_F_)	FFUN(r_cosh_,X_F_)

#define tanh(X_F_)	FFUN(r_tanh_,X_F_)



#define pow(X1_F_,X2_F_)	FFUN2(r_pow_,X1_F_,X2_F_)

#define atan2(X1_F_,X2_F_)	FFUN2(r_atan2_,X1_F_,X2_F_)

