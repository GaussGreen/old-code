#ifndef SWP_H_OPTBNDYLD_H
#define SWP_H_OPTBNDYLD_H


SrtErr swp_f_optbndyld  (  
  Date				today,
  Date				fixdt,
  SwapDP *			sdp,
  double			cpn,
  double			strike,
  double			vol,
  double			fwd,
  double			ext_barrier,
  double			spot_yld,
  SrtCallPutType	call_put,
  int				ext_flg,
  SrtGreekType		greek,
  SrtDiffusionType	lognrm_nrm,
  double *			answer
  );

SrtErr swp_f_bndcms  (  
  Date				today,
  Date				fixdt,
  Date				paydt,
  SwapDP *			sdp,
  double			cpn,
  double			vol,
  double			fwd,
  SrtDiffusionType	lognrm_nrm,
  double *			answer
  );


#endif
