
#include "SRT_H_ALL.H>
#include "pde_h_struct.h"

Err und_pde_lim(Date			event_date, 
				SrtUndPtr		und,
				SrtUndPtr		dom_und,
				double			max_time,
				long			min_node,
				long			min_num_mesh,
				SrtBasicPdeInf	*pde_inf);


Err srt_f_set_srvgs_pde_operators(Date				event_date,
							  String			dom_und_name,
							  String			und_name,
							  SrtBasicPdeInf	pde_inf,

							  double			strike,
							  SrtReceiverType	RecPay,
							  double			****diff_mat,
							  double			****conv_mat,
							  double			****div_tensor);
