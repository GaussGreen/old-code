#ifndef IMSL_KR_H

#define IMSL_KR_H



void	 *imsl_malloc();

void	 imsl_free();



void     imsl_e1psh();

void     imsl_e1pop();

void     imsl_e1csm();

void	 imsl_ercode();

void	 imsl_e1ucs();

void     imsl_e1usr();

void     imsl_e1pos();

void     imsl_e1mes();

void	 imsl_e1chk();

int      imsl_n1rty();

int      imsl_n1rcd();

int      imsl_n1rnof();



float    imsl_amach();

double   imsl_dmach();

void     imsl_umach();

char*    imsl_achar();

int      imsl_iachar();



	/* Chapter 1 */



void     imsl_lftrg();

void     imsl_dlftrg();

void     imsl_l2trg();

void     imsl_dl2trg();

void     imsl_lfsrg();

void     imsl_dlfsrg();



	/* Chapter 3 */



void     imsl_csint();

void     imsl_dcsint();

void     imsl_c2int();

void     imsl_dc2int();

void     imsl_c2dec();

void     imsl_dc2dec();

void     imsl_c1sor();

void     imsl_dc1sor();

float    imsl_csval();

double   imsl_dcsval();

float    imsl_ppder();

double   imsl_dppder();

void     imsl_p3der();

void     imsl_dp3der();

float    imsl_csder();

double   imsl_dcsder();

float    imsl_csitg();

double   imsl_dcsitg();

float    imsl_ppitg();

double   imsl_dppitg();

float    imsl_b2der();

double   imsl_db2der();

float    imsl_b3der();

double   imsl_db3der();

void     imsl_b4der();

void     imsl_db4der();

void     imsl_b2int();

void     imsl_db2int();

void     imsl_b3int();

void     imsl_db3int();

void     imsl_b4int();

void     imsl_db4int();

void     imsl_b5int();

void     imsl_db5int();

float    imsl_b2itg();

double   imsl_db2itg();

float    imsl_b3itg();

double   imsl_db3itg();

float    imsl_b4itg() ;

double   imsl_db4itg() ;

void     imsl_b5itg();

void     imsl_db5itg();

void     imsl_b2nak();

void     imsl_db2nak();

void     imsl_c1not();

void     imsl_dc1not();

void     imsl_b2opk();

void     imsl_db2opk();

void     imsl_b3opk() ;

void     imsl_db3opk() ;

void     imsl_b4opk();

void     imsl_db4opk();

void     imsl_crbrb();

void     imsl_dcrbrb();

void     imsl_l2lrb();

void     imsl_dl2lrb();

void     imsl_l2crb();

void     imsl_dl2crb();

void     imsl_l2trb();

void     imsl_dl2trb();

void     imsl_lfirb();

void     imsl_dlfirb();

void     imsl_lfsrb();

void     imsl_dlfsrb();

void     imsl_nr1rb();

void     imsl_dnr1rb();

void     imsl_stbsv();

void     imsl_dstbsv();

float    imsl_b22dr();

double   imsl_db22dr();

void     imsl_b32dr();

void     imsl_db32dr();

float    imsl_b22ig();

double   imsl_db22ig();

float    imsl_b32ig(); 

double   imsl_db32ig(); 

void     imsl_b22in();

void     imsl_db22in();

void     imsl_b32in();

void     imsl_db32in();

void     imsl_b42in();

void     imsl_db42in();

void     imsl_b2ls2(); 

void     imsl_db2ls2(); 

void     imsl_b3ls2();

void     imsl_db3ls2();

void     imsl_b4ls2();

void     imsl_db4ls2();

void     imsl_b5ls2();

void     imsl_db5ls2();

void     imsl_b6ls2();

void     imsl_db6ls2();

void     imsl_b7ls2();

void     imsl_db7ls2();

void     imsl_b2lsq();

void     imsl_db2lsq();

void     imsl_b3lsq();

void     imsl_db3lsq();

void     imsl_b4lsq();

void     imsl_db4lsq();

void     imsl_b5lsq();

void     imsl_db5lsq();

void     imsl_b6lsq();

void     imsl_db6lsq();

void     imsl_b2vls();

void     imsl_db2vls();

void     imsl_b3vls();

void     imsl_db3vls();

void     imsl_b4vls();

void     imsl_db4vls();

void     imsl_b5vls();

void     imsl_db5vls();

float    imsl_b6vls();

double   imsl_db6vls();

void     imsl_b7vls();

void     imsl_db7vls();

void     imsl_b8vls();

void     imsl_db8vls();

void     imsl_b2cpp();

void     imsl_db2cpp();

void     imsl_b3cpp();

void     imsl_db3cpp();

void     imsl_f2lsq();

void     imsl_df2lsq();

void     imsl_c2scv();

void     imsl_dc2scv();

void     imsl_c3scv(); 

void     imsl_dc3scv(); 

void     imsl_c4scv();

void     imsl_dc4scv();

void     imsl_c5scv();

void     imsl_dc5scv();

void     imsl_c2smh();

void     imsl_dc2smh();

void     imsl_c3smh();

void     imsl_dc3smh();

float    imsl_c4smh();

double   imsl_dc4smh();





	/* Chapter 4 */



void     imsl_qdng();

void     imsl_dqdng();

void     imsl_q3ng();

void     imsl_dq3ng();

void     imsl_q4ng();

void     imsl_dq4ng();



void     imsl_qdag();

void     imsl_dqdag();

void     imsl_q2ag();

void     imsl_dq2ag();

void     imsl_g2rcf();

void     imsl_dg2rcf();

float    imsl_g3rcf();

double   imsl_dg3rcf();

void     imsl_g4rcf();

void     imsl_dg4rcf();

void     imsl_g2rul();

void     imsl_dg2rul();

void     imsl_reccf();

void     imsl_dreccf();



	/* Chapter 5*/



void     imsl_ivprk();

void     imsl_divprk();

void     imsl_i2prk();

void     imsl_di2prk();

void     imsl_i3prk();

void     imsl_di3prk();



	/* Chapter 6 */



void     imsl_fftri();

void     imsl_fftrf();

void     imsl_fftrb();



	/* Level 1 BLAS */



int      imsl_isamax();

int      imsl_idamax();

#ifndef INLINE_BLAS

void     sscal();

void     dscal();

void     saxpy();

void     daxpy();

void     sswap();

void     dswap();

void     iset();

void     sset();

void     dset();

void     icopy();

void     scopy();

void     dcopy();

void     sadd();

void     dadd();

#endif

float    imsl_snrm2();

double   imsl_dnrm2();

float    imsl_sdot();

double   imsl_ddot();



	/* Level 9 BLAS */



void     imsl_sger();

void     imsl_dger();

void     imsl_strsv();

void     imsl_dtrsv();



	/* Chapter 9 */



void     trnrr_cf();

void     dtrnrr_cf();

void     imsl_crgrg();

void     imsl_dcrgrg();



	/* Chapter 10 */



void	 imsl_r2ivn();

void     imsl_c1wfr();



        /* Chapter 12 */



void	 imsl_rnses();

void	 imsl_r1int();

void	 imsl_r1clk();



	/* Utilities*/



void     svrgp_c();

void     dsvrgp_c();

void     svrgn_c();

void     dsvrgn_c();

void     svign_c();



	/* SFUN/LIBRARY*/



int      imsl_inits();

int      imsl_initds();

float    imsl_csevl();

double   imsl_dcsevl();



float    imsl_erfce();

double   imsl_derfce();



float    imsl_beta();

float    imsl_betai();

float    imsl_gami();

float    imsl_gamit();

float    imsl_gamma();

float    imsl_gamr();



float    imsl_bsi0e();

float    imsl_bsi1e();

float    imsl_bsk0e();

float    imsl_bsk1e();

double	 imsl_dbsi0e();

double	 imsl_dbsi1e();

double	 imsl_dbsk0e();

double	 imsl_dbsk1e();



float    imsl_anorin();

float    imsl_anordf();

float    imsl_chidf();

float    imsl_chiin();

float    imsl_fdf();

float    imsl_fin();

float    imsl_tdf();

float    imsl_tin();

float    imsl_gamdf();

float    imsl_bindf();

float    imsl_hypdf();

float    imsl_poidf();

float    imsl_betin();

float    imsl_betdf();



float    imsl_alnrel();

float    imsl_albeta();

float    imsl_alngam();

void     imsl_algams();



void     imsl_r9gaml();

float    imsl_r9gmit();

float    imsl_r9lgic();

float    imsl_r9lgit();

float    imsl_gamma();

double   imsl_dgamma();

float    imsl_r9lgmc();

double   imsl_d9lgmc();



float    imsl_carg();

double   imsl_zarg();

f_complex imsl_clngam();

f_complex imsl_clnrel();

f_complex imsl_c9lgmc();



	/* AT Tests */



void     imsl_ststn();

void	 imsl_rititm();

void     imsl_pasckd();

void     imsl_pascks();

void     imsl_pasckt();

void     imsl_abnrm();

void     imsl_diffs();

void     imsl_wrtrl();

void     imsl_e1ucs();



       /*  Stat error checking routines */

void            imsl_c12ile();

void            imsl_c1iarg();



      /*  Transpose  */

void            imsl_i_m1ran();

void            imsl_f_m1ran();

void            imsl_c_m1ran();

void            imsl_prime();



     /*  The writing routines  */

void            imsl_write_controller();

Mint            imsl_write_initialize();

void            imsl_write_labels();

void            imsl_write_title();

Mchar          *imsl_write_conversion();

void            imsl_c1nter();

void            imsl_w5rrl_f();

void            imsl_w6rrl();

Mchar          *imsl_w7rrl();

void            imsl_w8rrl();

void            imsl_w12rl();

void            imsl_wropt();

void            imsl_w1opt();

void            imsl_write_line();

Mchar          *imsl_w1iss();

Mchar          *imsl_fmtx();

void            imsl_write_format();



    /*  Stat bla (checks for missing) */

void            imsl_i1max();



#endif /* IMSL_KR_H */



