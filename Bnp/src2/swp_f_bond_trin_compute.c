/************************************************************************/
/* MODULE NAME : 	BOND TRINOMIAL					*/
/* Version 1.00  	4/11/1994					*/
/* Author		A.P.  , O.VE.  , A.B.  , HdL */
/************************************************************************/
/********
   Last modified: Dec 2  , 1994
   Atuhor : K L Chau
   Reason : use new Srtcrv structure
*********/

#include "opfnctns.h"
#include "srt_h_all.h"

#define DEBUG_VEGA 0
#define DEBUG_PRINT 0
#define DEBUG_GREEKS 0
#define H_DEBUG 0
#define H_DEBUG_1 0
#define YIELD_DEBUG 0
#define BUG_DEBUG 0

/**********************functions declaration***********************/

double fwd_clean_fct(Ddate dstart, Ddate fut, String swap, String repo,
                     Deal contract);
void tree_init(Deal contract, Model *model, Step *tree_steps, double *vega);
void tree_choose(Deal *contract, Model *model, Func *func);
void calc_proba(Model model, Step *step);
void shift_all_following_steps(Model model, int index, Step *tree_steps,
                               int *step_num);
void add_step(Deal contract, Model model, int index_tree, int index_date,
              Step *tree_steps);
void date_treatment(Deal contract, Model model, Model model_fx,
                    Step *tree_steps, Step *tree_steps_fx, int *step_num,
                    int *step_fx_num);
void init_step(int n, Deal contract, Model model, Model model_fx, Step *step,
               Step *step_fx);
double q2(double alpha, double a, double b);
void smoothing(Deal contract, Model model, Step *step, Node *node);
double endtree_european(Deal contract, Model model, Step step, Node node);
double endtree_partial(Deal contract, Model model, Step step, Node node);
double payoff_european(Deal contract, Model model, Step step, Node node);
double payoff_american(Deal contract, Model model, Step step, Node node);
double endtree_shout_lock(Deal contract, Model model, Step step, Node node);
double payoff_shout_lock(Deal contract, Model model, Step step, Node node);
double endtree_shout_reset(Deal contract, Model model, Step step, Node node);
double payoff_shout_reset(Deal contract, Model model, Step step, Node node);
double endtree_both_shout(Deal contract, Model model, Step step, Node node);
double payoff_both_shout(Deal contract, Model model, Step step, Node node);
double payoff_instal(Deal contract, Model model, Step step, Node node);
double payoff_extinguish(Deal contract, Model model, Step step, Node node);
double endtree_extinguish(Deal contract, Model model, Step step, Node node);
double payoff_mid_atlantic(Deal contract, Model model, Step step, Node node);
double endtree_you_choose(Deal contract, Model model, Step step, Node node);
double endtree_digit(Deal contract, Model model, Step step, Node node);
double payoff_amer_digit(Deal contract, Model model, Step step, Node node);
double endtree_pay_digit(Deal contract, Model model, Step step, Node node);
double payoff_pay_digit(Deal contract, Model model, Step step, Node node);

/*** function to calculate the forward price ***/
/*** at a date which can be a double	     ***/

double fwd_clean_fct(Ddate dstart, Ddate fut, String swap, String repo,
                     Deal contract) {
  int nc, i;
  double dirty_price, fwd_clean_price, acc, fwd_acc_int;
  double sum_coup = 0.0, c, df, theta;
  DateList list;
  SwapDP pf;

  pf = contract.p;

  fut += (double)contract.d_set;
  /****** Here   , we take into account the nb of day settlement for the
          bond in all the calculations of the fwd price	********/

  acc = acc_int_fct(contract.p, contract.first_coupon);
  dirty_price = contract.spot + acc;

  list = SwapDP_to_DateList(&pf, NO_BUSDAY_CONVENTION);
  c = contract.first_coupon;
  nc = 0;
  while (list.date[nc + 1] <= fut)
    nc++;
  for (i = 1; i <= nc; i++) {
    sum_coup += (c / contract.p.compd) /
                df_crv((Ddate)(list.date[i]), (Ddate)fut, lookup_curve(swap));
    c = contract.coupon;
  }

  pf.start = (long)fut;
  theta = (fut - (double)pf.start) / 365.0;
  fwd_acc_int = acc_int_fct(pf, c);
  fwd_acc_int += theta * contract.coupon;
  df = df_crv((Ddate)dstart, (Ddate)fut, lookup_curve(repo), 1);
  fwd_clean_price = dirty_price / df - sum_coup - fwd_acc_int;
  srt_free(list.date);
  return (fwd_clean_price);
}

/*************** choose the type of option *******************/

void tree_choose(Deal *contract, Model *model, Func *func) {

  if (contract->type == EUROPEAN) {
    func->PAY_OFF = payoff_european;
    func->END_TREE = endtree_european;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == AMERICAN) {
    func->PAY_OFF = payoff_american;
    func->END_TREE = endtree_european;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == SHOUT_LOCK) {
    func->PAY_OFF = payoff_shout_lock;
    func->END_TREE = endtree_shout_lock;
    model->h = (double)((contract->sub_period - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == SHOUT_RESET) {
    func->PAY_OFF = payoff_shout_reset;
    func->END_TREE = endtree_shout_reset;
    model->h = (double)((contract->sub_period - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == BOTH_SHOUT) {
    func->PAY_OFF = payoff_both_shout;
    func->END_TREE = endtree_both_shout;
    model->h = (double)((contract->sub_period - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == INSTALMENT) {
    func->PAY_OFF = payoff_instal;
    func->END_TREE = endtree_european;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
  } else if (contract->type == EXTINGUISH) {
    func->PAY_OFF = payoff_extinguish;
    func->END_TREE = endtree_extinguish;
    model->h = (double)((contract->sub_period - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == MID_ATLANTIC) {
    func->PAY_OFF = payoff_mid_atlantic;
    func->END_TREE = endtree_european;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
  } else if (contract->type == YOU_CHOOSE) {
    /***	contract->sub_period will be the choose date **/
    /***             ->strike                 call strike **/
    /***             ->barrier                put  strike **/

    func->PAY_OFF = payoff_european;
    func->END_TREE = endtree_you_choose;
    model->h = (double)((contract->sub_period - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == EURO_DIGIT) {
    func->PAY_OFF = payoff_european;
    func->END_TREE = endtree_digit;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == AMER_DIGIT) {
    /****	CAUTION			****/
    /**** up gives an up and out   	****/

    func->PAY_OFF = payoff_amer_digit;
    func->END_TREE = endtree_digit;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  } else if (contract->type == PAY_DIGIT) {
    /****	CAUTION			****/
    /**** up gives an up and in 	****/

    func->PAY_OFF = payoff_pay_digit;
    func->END_TREE = endtree_pay_digit;
    contract->sub_period = contract->opt_mat;
    model->h = (double)((contract->opt_mat - contract->today) /
                        (model->step_num * 365.0));
    contract->size_date_cpn = 0;
  }
}

/************* treatment for instalment dates *****************/
/***** add steps  , find position of instalment dates ***********/

void shift_all_following_steps(Model model, int index, Step *step,
                               int *step_num) {
  int i;

  for (i = *step_num; i >= index; i--) {
    step[i + 1] = step[i];
  }
  *step_num += 1;
}

void add_step(Deal contract, Model model, int index_tree, int index_date,
              Step *step) {
  double t_n_y;
  SrtCrvPtr crv;

  crv = lookup_curve(model.swap);

  step[index_tree].df_swap =
      df_crv((Ddate)(contract.today), (Ddate)(contract.date[index_date]), crv);
  t_n_y = (double)((contract.date[index_date] - contract.today) / 365.0);
  /*	useless  = interp(	&model.date_\cumvol[0][0]  ,
                                  &model.date_cumvol[1][0]  ,
                                  model.size_tsv  ,
                                  (double)(t_n_y)  ,
                                  0  ,
                                  &temp);
  */
  step[index_tree].cum_vol = eq_cum_vol_func((double)t_n_y, SRT_YES, model.ts);
  step[index_tree].time_start = t_n_y;
  step[index_tree].flag_date = 1;
  step[index_tree].cpn = contract.cpn[index_date];
  step[index_tree].expected_spot = fwd_clean_fct(
      (double)contract.p.start, (double)contract.today + t_n_y * 365.0,
      model.swap, model.repo, contract);
  DEBUG_PRINT &&printf("add step:%d  exp_spot:%f \n", index_tree,
                       step[index_tree].expected_spot);

  /******** wil be called again for quantos *******/
}

void date_treatment(Deal contract, Model model, Model model_fx, Step *step,
                    Step *step_fx, int *step_num, int *step_fx_num) {
  int i, j;

  if ((contract.type == INSTALMENT) || (contract.type == MID_ATLANTIC)) {
    for (i = 0; i < contract.size_date_cpn; i++) {
      for (j = 0; j <= *step_num; j++) {
        if (contract.date[i] <=
            (double)(365.0 * step[j].time_start + contract.today))
          break;
      }

      if ((j > 0) && ((j < *step_num) &&
                      (contract.date[i] < (double)(365.0 * step[j].time_start +
                                                   contract.today)))) {
        shift_all_following_steps(model, j, step, step_num);
        add_step(contract, model, j, i, step);

        if (contract.quanto == 1) {
          shift_all_following_steps(model_fx, j, step_fx, step_fx_num);
          add_step(contract, model_fx, j, i, step_fx);
          step[j].expected_spot *= exp(
              contract.corr_fx * sqrt(step[j].cum_vol * step_fx[j].cum_vol));
          step[j].df_swap = step_fx[j].df_swap;
        }

      } else {
        if ((j > 0) &&
            ((j < *step_num) &&
             (contract.date[i] ==
              (double)(365.0 * step[j].time_start + contract.today)))) {
          step[j].flag_date = 1;
          step[j].cpn = contract.cpn[i];
          DEBUG_PRINT &&printf("date_cpn on step:%d \n", j);
        }
      }
    }
  }
}

/****************** initialization ****************************/

void tree_init(Deal contract, Model *model, Step *step, double *vega) {
  int i, j;
  double *vol_inst, t_n_end;
  double t_n, t_n_y, vol_max, interm;

  /***************** memory allocation ****************/
  model->date_cumvol = dmatrix(0, 1, 0, model->size_tsv - 1);
  vol_inst = dvector(0, model->size_tsv - 1);

  /****** determination of delta_lns + calculation of cumvol *****/

  model->date_cumvol[0][0] =
      (double)((model->tsv[0][0] - contract.today) / 365.0);
  model->date_cumvol[1][0] =
      model->tsv[1][0] * model->tsv[1][0] * model->date_cumvol[0][0];
  vol_max = model->tsv[1][0];

  for (i = 1; i < model->size_tsv; i++) {
    model->date_cumvol[0][i] =
        (double)((model->tsv[0][i] - contract.today) / 365);
    model->date_cumvol[1][i] =
        model->tsv[1][i] * model->tsv[1][i] * model->date_cumvol[0][i];
    vol_inst[i] =
        sqrt((model->date_cumvol[1][i] - model->date_cumvol[1][i - 1]) /
             (model->date_cumvol[0][i] - model->date_cumvol[0][i - 1]));

    if (vol_inst[i] > vol_max)
      vol_max = vol_inst[i];
  }

  for (i = 1; i < model->size_tsv; i++) {
    model->mat_vol_nb = i;
    if (model->date_cumvol[0][i] > (double)contract.opt_mat)
      break;
  }

  model->delta_lns = sqrt(3 * model->h) * vol_max;
  interm = model->delta_lns;

  DEBUG_VEGA &&printf("1 vega=%f \n", *vega);
  if ((contract.ret_type == 3) || (contract.ret_type == 4)) {
    if (*vega > EPSILON)
      model->delta_lns = *vega;
    else
      *vega = interm;
    DEBUG_VEGA &&printf("2 vega=%f \n", *vega);
  }

  model->diff_up = contract.spot * (exp(model->delta_lns) - 1);
  model->diff_do = contract.spot * (1 - exp(-model->delta_lns));
  DEBUG_VEGA &&printf("DIFF_UP=%f \n", model->diff_up);

  /*************** initialization of df and vol ***************/

  t_n_end = (double)((contract.opt_mat - contract.today) / 365.0);
  /*	useless_end  = interp(	&model->date_cumvol[0][0]  ,
                                  &model->date_cumvol[1][0]  ,
                                  model->size_tsv  ,
                                  (double)(t_n_end)  ,
                                  0  ,
                                  &temp_end);
  */
  model->cum_vol_opt = eq_cum_vol_func((double)t_n_end, SRT_YES, model->ts);
  model->df_swap_all =
      df_crv((Ddate)(contract.today), (Ddate)(contract.opt_mat),
             lookup_curve(model->swap));

  t_n = (double)(contract.today);

  for (j = 0; j <= model->step_num; j++) /***  at each step  ***/
  {
    t_n_y = (double)((t_n - contract.today) / 365.0);
    /*        	useless_end  = interp(	&model->date_cumvol[0][0]  ,
                                            &model->date_cumvol[1][0]  ,
                                            model->size_tsv  ,
                                            (double)(t_n_y)  ,
                                            0  ,
                                            &temp_end);
    */
    step[j].cum_vol = eq_cum_vol_func((double)t_n_y, SRT_YES, model->ts);
    step[j].time_start = t_n_y;
    step[j].flag_date = 0;
    step[j].df_swap =
        df_crv((Ddate)(contract.today), (Ddate)t_n, lookup_curve(model->swap));

    t_n += model->h * 365.0;
  }

  free_dvector(vol_inst, 0, model->size_tsv - 1);

  /*******  for extinguish option hit_level position  *******/
  if ((contract.type == EXTINGUISH) || (contract.type == AMER_DIGIT) ||
      (contract.type == EURO_DIGIT) || (contract.type == PAY_DIGIT))

  {

    model->hit_level =
        (int)((2 * contract.do_up - 1) *
              floor(log(contract.barrier / contract.spot) *
                    (2.0 * contract.do_up - 1.0) / model->delta_lns));

    /* Hit level in case of yield barriers  */
    /* will be defined in step_init (later) */

    DEBUG_PRINT &&printf("barrier=%f  spot=%f  ", contract.barrier,
                         contract.spot);
    DEBUG_PRINT &&printf("delta=%f  \n", model->delta_lns);
    H_DEBUG &&printf("hit_level defined : %d \n", model->hit_level);
  }
}

/******************** calculation of probabilities *****************/

void calc_proba(Model model, Step *step) {
  double m, var;
  Step *test = step;
  test++;

  m = log(test->expected_spot / step->expected_spot);
  m -= step->vol_h / 2;
  step->drift_h = m;
  var = step->vol_h;
  step->d_index = (int)floor(m / model.delta_lns + 0.5);
  m -= step->d_index * model.delta_lns;

  step->p_d = (m * m + var - model.delta_lns * m) /
              (2 * model.delta_lns * model.delta_lns);
  step->p_u = (m * m + var + model.delta_lns * m) /
              (2 * model.delta_lns * model.delta_lns);
  step->p_m = 1 - step->p_u - step->p_d;

  DEBUG_PRINT &&printf("proba up:%f  mid:%f  down:%f \n", step->p_u, step->p_m,
                       step->p_d);
}

/**********************************************************/
/***********   initialization of step n   *****************/
/**********************************************************/

void init_step(int n, Deal contract, Model model, Model model_fx, Step *step,
               Step *step_fx) {
  SwapDP p_temp;

  step[n].df_swap_end = model.df_swap_all / step[n].df_swap;
  step[n].time_end = (double)((contract.opt_mat - contract.today) / 365.0 -
                              step[n].time_start);
  step[n].expected_spot =
      fwd_clean_fct((double)contract.p.start,
                    (double)contract.today + step[n].time_start * 365.0,
                    model.swap, model.repo, contract);

  if (contract.yield_strike_flag == 1) {
    p_temp = contract.p;
    p_temp.start = (int)(floor(contract.p.start + step[n].time_start * 365));
    step[n].clean_strike = 100 * clean_price_fct(p_temp, contract.coupon / 100,
                                                 contract.yield_strike / 100, 1,
                                                 contract.coupon / 100);
  }

  if (contract.yield_barrier_flag == 1) {
    p_temp = contract.p;
    p_temp.start = (int)(floor(contract.p.start + step[n].time_start * 365));
    step[n].clean_barrier = 100 * clean_price_fct(p_temp, contract.coupon / 100,
                                                  contract.yield_barrier / 100,
                                                  1, contract.coupon / 100);

    if ((contract.type == EXTINGUISH) || (contract.type == AMER_DIGIT) ||
        (contract.type == EURO_DIGIT) || (contract.type == PAY_DIGIT)) {
      step[n].yield_hit_level =
          (int)((2 * contract.do_up - 1) *
                floor(log(step[n].clean_barrier / contract.spot) *
                      (2.0 * contract.do_up - 1.0) / model.delta_lns));
    }

  } /*** end of yield barriers **/

  YIELD_DEBUG &&printf("clean strike=%f barrier=%f \n", step[n].clean_strike,
                       step[n].clean_barrier);
  YIELD_DEBUG &&printf("time start =%f \n", step[n].time_start);

  if (contract.quanto == 1) {
    step[n].df_swap_end = model_fx.df_swap_all / step_fx[n].df_swap;
    step[n].expected_spot *=
        exp(contract.corr_fx * sqrt(step[n].cum_vol * step_fx[n].cum_vol));

    H_DEBUG &&printf("init steps spot expected old =%f \n",
                     step[n].expected_spot);
    H_DEBUG &&printf("corr =%f cum_vol =%f cum_fx =%f \n ", contract.corr_fx,
                     step[n].cum_vol, step_fx[n].cum_vol);
    H_DEBUG &&printf("init steps spot  new =%f \n", step[n].expected_spot);
  }

  if (n < model.step_num) {
    step[n].vol_h = step[n + 1].cum_vol - step[n].cum_vol;
    step[n].vol =
        sqrt((model.cum_vol_opt - step[n].cum_vol) / step[n].time_end);

    if (contract.quanto == 1) { /*** uses the fx swap **/
      step[n].df_swap_h = step_fx[n + 1].df_swap / step_fx[n].df_swap;
    } else {
      step[n].df_swap_h = step[n + 1].df_swap / step[n].df_swap;
    }

    calc_proba(model, &step[n]);
    /**** OK for B.S.+ where probas are equal at each node	****/
  }

  if ((n == model.step_num) &&
      (contract.sub_period !=
       contract.opt_mat)) { /****	Still need a vol for B-S at the end of the
                               tree	****/
    step[n].vol =
        sqrt((model.cum_vol_opt - step[n].cum_vol) / step[n].time_end);
  }
}

/*************** smoothing for Barrier option : change var***************/

double q2(double alpha, double a, double b) {
  double result;

  result = exp(a * alpha + 0.5 * b * b * alpha * alpha);
  result = result * norm((a + alpha * b * b) / b);
  return (result);
}

void smoothing(Deal contract, Model model, Step *step, Node *node) {
  double u, d, do_up;
  double m, m1, m2, ft, fe, absorb_prob, mu_h;
  double absorb_moment[1], absorb_moments[1];
  Step *test;

  if (((contract.type == EXTINGUISH) || (contract.type == AMER_DIGIT) ||
       (contract.type == EURO_DIGIT) || (contract.type == PAY_DIGIT)) &&
      (((node->position == model.hit_level) &&
        (contract.yield_barrier_flag == 0)) ||
       ((node->position == step->yield_hit_level) &&
        (contract.yield_barrier_flag == 1))))

  {
    if (contract.yield_barrier_flag == 1)
      contract.barrier = step->clean_barrier;

    do_up = ((contract.do_up == 0) ? -1 : 1);
    u = exp(model.delta_lns);
    d = 1 / u;
    test = step;
    test++;
    m = log(contract.barrier / node->asset);
    mu_h = log(test->expected_spot / step->expected_spot);
    mu_h -= step->vol_h / 2;
    m1 = (m - mu_h) * do_up;
    m2 = -(m + mu_h) * do_up;
    ft = sqrt(step->vol_h);
    fe = exp(2 * mu_h * m / step->vol_h);

    absorb_prob = norm(m1 / ft) - fe * norm(m2 / ft);

    absorb_moment[0] =
        contract.barrier * (q2(-do_up, m1, ft) - fe * q2(-do_up, m2, ft));
    absorb_moment[0] /= absorb_prob;
    absorb_moments[0] =
        absorb_moment[0] * absorb_prob + contract.barrier * (1 - absorb_prob);
    absorb_moments[0] /= node->asset;

    absorb_moment[1] = contract.barrier * contract.barrier *
                       (q2(-2 * do_up, m1, ft) - fe * q2(-2 * do_up, m2, ft)) /
                       absorb_prob;

    absorb_moments[1] = absorb_moment[1] * absorb_prob +
                        contract.barrier * contract.barrier * (1 - absorb_prob);
    absorb_moments[1] /= node->asset * node->asset;

    if (do_up == 1) /* Up */
    {
      node->p_d = (absorb_moments[1] - absorb_moments[0] * absorb_moments[0]) +
                  (contract.barrier / node->asset - absorb_moments[0]) *
                      (1 - absorb_moments[0]);
      node->p_d /= (1 - d) * (contract.barrier / node->asset - d);

      node->p_m = (absorb_moments[1] - absorb_moments[0] * absorb_moments[0]) +
                  (contract.barrier / node->asset - absorb_moments[0]) *
                      (d - absorb_moments[0]);
      node->p_m /= (d - 1) * (contract.barrier / node->asset - 1);

      node->p_u = 1 - node->p_m - node->p_d;
    } else {
      node->p_d = (absorb_moments[1] - absorb_moments[0] * absorb_moments[0]) +
                  (u - absorb_moments[0]) * (1 - absorb_moments[0]);
      node->p_d /= (1 - contract.barrier / node->asset) *
                   (u - contract.barrier / node->asset);

      node->p_m = (absorb_moments[1] - absorb_moments[0] * absorb_moments[0]) +
                  (u - absorb_moments[0]) *
                      (contract.barrier / node->asset - absorb_moments[0]);
      node->p_m /= (contract.barrier / node->asset - 1) * (u - 1);

      node->p_u = 1 - node->p_m - node->p_d;
    }

    DEBUG_PRINT &&printf("niveau HIT ! node->asset=%f\n", node->asset);
    DEBUG_PRINT &&printf("proba up:%f  mid:%f  down:%f \n", node->p_u,
                         node->p_m, node->p_d);
  }
}

/****************** the TREE **************************/

double bond_trinomial_fct(Deal contract, Model model, Model model_fx,
                          double *vega) {

  int i, k, n, d_index, d0;
  double price, delta, gamma, theta_1d, down;
  double p_up, p_do, p_mi, disc_fact;
  double *prem_prev;
  double *prem_cur;
  Node node;
  Func func;
  Step *steps;
  Step *steps_fx;
  SwapDP pf, pfwd;
  int opt_num;
  double temp;

  BUG_DEBUG &&printf("1 \n");

  /**************** init ********************/

  prem_prev = calloc(2 * model.step_num + 50, sizeof(double));
  if (!prem_prev)
    printf("erreur chez prem_prev\n");

  prem_cur = calloc(2 * model.step_num + 50, sizeof(double));
  if (!prem_cur)
    printf("erreur chez prem_cur\n");

  steps = calloc(model.step_num + 50, sizeof(Step));
  if (!steps)
    printf("erreur chez steps\n");

  tree_choose(&contract, &model, &func);

  tree_init(contract, &model, steps, vega);
  DEBUG_VEGA &&printf("3 vega=%f \n", *vega);

  /*******************************************************/
  /** optimization of the number of steps for 1 barrier **/
  /******** 	only for euro digits 		********/
  /*******************************************************/

  if (contract.type == EURO_DIGIT) {
    temp = fabs(log(contract.barrier / contract.spot));
    temp /= model.delta_lns * sqrt(model.step_num);

    H_DEBUG_1 &&printf(" temp  =%f  \n", temp);

    for (i = 1; i < model.step_num; i++) {
      opt_num = (int)floor(i * i * temp * temp) + contract.call_put;

      H_DEBUG_1 &&printf(" opt_num =%d \n", opt_num);

      if (opt_num > model.step_num) {
        model.step_num = opt_num;
        break;
      }
    }

    H_DEBUG_1 &&printf(" new number of steps  =%d \n", model.step_num);

    /****	FREE allocated memory before reallocation with new step nb
     * ****/

    srt_free(prem_cur);
    srt_free(prem_prev);
    srt_free(steps);
    free_dmatrix(model.date_cumvol, 0, 1, 0, model.size_tsv - 1);

    prem_prev = calloc(2 * model.step_num + 50, sizeof(double));
    if (!prem_prev)
      printf("erreur chez prem_prev\n");

    prem_cur = calloc(2 * model.step_num + 50, sizeof(double));
    if (!prem_cur)
      printf("erreur chez prem_cur\n");

    steps = calloc(model.step_num + 50, sizeof(Step));
    if (!steps)
      printf("erreur chez steps\n");

    tree_choose(&contract, &model, &func);
    tree_init(contract, &model, steps, vega);

  } /**** end of if EURO_DIGIT ********/

  model.fwd_opt_mat = fwd_clean_fct((double)contract.p.start, contract.opt_mat,
                                    model.swap, model.repo, contract);

  /**** for quantos ****/

  if (contract.quanto == 1) {
    model_fx.step_num = model.step_num;
    model_fx.h = model.h;

    H_DEBUG &&printf(" \n starts quantos with corr =%f \n ", contract.corr_fx);

    steps_fx = calloc(model_fx.step_num + 50, sizeof(Step));
    if (!steps_fx)
      printf("erreur chez steps_fx\n");

    tree_init(contract, &model_fx, steps_fx, vega);

    H_DEBUG &&printf("old forward mat =%f \n", model.fwd_opt_mat);

    model.fwd_opt_mat *= exp(
        contract.corr_fx * sqrt(model.date_cumvol[1][model.mat_vol_nb] *
                                model_fx.date_cumvol[1][model_fx.mat_vol_nb]));

    H_DEBUG &&printf("    new forward mat =%f \n", model.fwd_opt_mat);
  }

  date_treatment(contract, model, model_fx, steps, steps_fx, &model.step_num,
                 &model_fx.step_num);

  pf = contract.p;
  pfwd = pf;

  /**************** extremity of the tree ********************/

  init_step(model.step_num, contract, model, model_fx, steps, steps_fx);

  node.asset = contract.spot * exp(model.delta_lns * (model.step_num + 1));
  down = exp(-model.delta_lns);

  /**** everything is already changed for quantos ****/

  node.forward =
      node.asset * model.fwd_opt_mat / steps[model.step_num].expected_spot;

  /**** in case of strike or barrier in yield ****/
  if (contract.yield_strike_flag == 1)
    contract.strike = steps[model.step_num].clean_strike;
  if (contract.yield_barrier_flag == 1)
    contract.barrier = steps[model.step_num].clean_barrier;

  for (i = 0; i <= 2 * model.step_num + 2; i++) {
    node.position = model.step_num + 1 - i;

    /*****	the node position lies between -(n+1) and +(n+1)	*****/
    /*****	whereas in the bitrin it is   0 to 2*(n+1)		*****/

    node.intrinsic =
        (node.asset - contract.strike) * (1 - 2 * contract.call_put);
    prem_prev[i] = func.END_TREE(contract, model, steps[model.step_num], node);
    node.asset *= down;
    node.forward *= down;
  }

  /**************** beginning the backwardization ****************/

  /**********************/
  for (n = model.step_num; n >= 1; n--) /***  for each step  **/
  {                                     /**********************/

    init_step(n - 1, contract, model, model_fx, steps, steps_fx);

    /****	initialize the first node at the top of the tree	****/
    node.asset = contract.spot * exp(model.delta_lns * n);
    node.position = n;

    /**** in case of strike or barrier in yield ****/
    if (contract.yield_strike_flag == 1)
      contract.strike = steps[model.step_num].clean_strike;
    if (contract.yield_barrier_flag == 1)
      contract.barrier = steps[model.step_num].clean_barrier;

    /***** 	top node : +n  *****/
    d_index = steps[n - 1].d_index;
    d0 = d_index; /* when  drift > up   , for node +n   , d0>0   */
    if (d0 >= 0)
      d0 = 0; /* leads to undefined node because node -1 */
              /* not defined - NEVER SUPPOSE CASE D0>1   */
    p_up = steps[n - 1].p_u;
    p_do = steps[n - 1].p_d;
    p_mi = steps[n - 1].p_m;
    disc_fact = steps[n - 1].df_swap_h;

    DEBUG_PRINT &&printf("node=%d ", node.position);

    if ((((contract.type == EXTINGUISH) || (contract.type == PAY_DIGIT) ||
          (contract.type == AMER_DIGIT)) ||
         ((contract.type == EURO_DIGIT) && (n == model.step_num))) &&
        (((node.position == model.hit_level) &&
          (contract.yield_barrier_flag == 0)) ||
         ((node.position == steps[n - 1].yield_hit_level) &&
          (contract.yield_barrier_flag == 1)))) {

      smoothing(contract, model, steps + n - 1, &node);

      node.disc_premium = node.p_u * prem_prev[0 - d0];
      node.disc_premium += node.p_m * prem_prev[1 - d0];
      node.disc_premium += node.p_d * prem_prev[2 - d0];
      node.disc_premium *= disc_fact;

      H_DEBUG &&printf("smooths with =%f =%f =%f \n", node.p_u, node.p_m,
                       node.p_d);
      H_DEBUG &&printf(" old smooth with =%f =%f =%f \n", p_up, p_mi, p_do);

    } else {
      node.disc_premium = p_up * prem_prev[0 - d0];
      node.disc_premium += p_mi * prem_prev[1 - d0];
      node.disc_premium += p_do * prem_prev[2 - d0];
      node.disc_premium *= disc_fact;
    }

    node.forward = node.asset * model.fwd_opt_mat / steps[n - 1].expected_spot;
    node.intrinsic =
        (node.asset - contract.strike) * (1 - 2 * contract.call_put);

    prem_cur[0] = func.PAY_OFF(contract, model, steps[n - 1], node);

    /*****  end of top node   *****/

    node.asset *= down;
    node.forward *= down;

    /**********************************/
    for (i = 1; i <= 2 * n - 1; i++) /*** for each node of the step  ***/
    {                                /*** except extremes for d0 pb  ***/
                                     /**********************************/
      node.position--;

      DEBUG_PRINT &&printf("%d ", node.position);
      H_DEBUG_1 &&printf(" position =%d \n", node.position);

      smoothing(contract, model, steps + n - 1, &node);

      if ((((contract.type == EXTINGUISH) || (contract.type == PAY_DIGIT) ||
            (contract.type == AMER_DIGIT)) ||
           ((contract.type == EURO_DIGIT) && (n == model.step_num) &&
            (fabs(node.asset - contract.barrier) > 0.00001))) &&
          (((node.position == model.hit_level) &&
            (contract.yield_barrier_flag == 0)) ||
           ((node.position == steps[n - 1].yield_hit_level) &&
            (contract.yield_barrier_flag == 1))))

      {
        node.disc_premium = node.p_u * prem_prev[i - d_index];
        node.disc_premium += node.p_m * prem_prev[i + 1 - d_index];
        node.disc_premium += node.p_d * prem_prev[i + 2 - d_index];
        node.disc_premium *= disc_fact;

        H_DEBUG_1 &&printf(" smoothed probas up =%f mid  =%f\n", node.p_u,
                           node.p_m);
        H_DEBUG_1 &&printf(" smoothed probas prem =%f \n", node.disc_premium);

      } else {
        node.disc_premium = p_up * prem_prev[i - d_index];
        node.disc_premium += p_mi * prem_prev[i + 1 - d_index];
        node.disc_premium += p_do * prem_prev[i + 2 - d_index];
        node.disc_premium *= disc_fact;

        H_DEBUG_1 &&printf("  probas up =%f mid  =%f\n", p_up, p_mi);
        H_DEBUG_1 &&printf("  prem =%f \n", node.disc_premium);
      }

      node.intrinsic =
          (node.asset - contract.strike) * (1 - 2 * contract.call_put);

      prem_cur[i] = func.PAY_OFF(contract, model, steps[n - 1], node);

      node.asset *= down;
      node.forward *= down;
    }

    node.position--;

    /****** down node : -n ******/

    DEBUG_PRINT &&printf("-n=%d   \n", node.position);

    smoothing(contract, model, steps + n - 1, &node);

    if (d_index > 0.0)
      d0 = d_index;
    else
      d0 = 0; /**** do not go to undefined node	*****/

    if (((contract.type == EXTINGUISH) || (contract.type == PAY_DIGIT) ||
         (contract.type == AMER_DIGIT) ||
         ((contract.type == EURO_DIGIT) && (n == model.step_num))) &&
        (((node.position == model.hit_level) &&
          (contract.yield_barrier_flag == 0)) ||
         ((node.position == steps[n - 1].yield_hit_level) &&
          (contract.yield_barrier_flag == 1))))

    {
      node.disc_premium = node.p_u * prem_prev[2 * n - d0];
      node.disc_premium += node.p_m * prem_prev[2 * n + 1 - d0];
      node.disc_premium += node.p_d * prem_prev[2 * n + 2 - d0];
      node.disc_premium *= disc_fact;
    } else {
      node.disc_premium = p_up * prem_prev[2 * n - d0];
      node.disc_premium += p_mi * prem_prev[2 * n + 1 - d0];
      node.disc_premium += p_do * prem_prev[2 * n + 2 - d0];
      node.disc_premium *= disc_fact;
    }

    node.intrinsic =
        (node.asset - contract.strike) * (1 - 2 * contract.call_put);

    prem_cur[2 * n] = func.PAY_OFF(contract, model, steps[n - 1], node);

    /****** end of node -n ******/

    /**************************************************/
    /***	CAUTION !!!				***/
    /***	node varies between +n to -n whereas	***/
    /***	prem_cur is indexed from 0 to 2*n   	***/
    /**************************************************/

    /***	prem_prev not changed for n=1 because of greeks calculation
     * ***/

    if (n > 1) {
      for (k = 0; k <= 2 * n; k++)
        prem_prev[k] = prem_cur[k];
    }

    if (n == 1)
      price = prem_cur[1];
  }

  /****	Determination of delta and gamma in the tree	****/
  /**** 	using the following calculation is a pb because ****/
  /****	delta_lns is too big only for barrier options	****/

  delta = (prem_cur[0] - prem_cur[2]) / (model.diff_up + model.diff_do);
  gamma = (prem_cur[0] - prem_cur[1]) / model.diff_up;
  gamma -= (prem_cur[1] - prem_cur[2]) / model.diff_do;
  gamma /= 0.5 * (model.diff_up + model.diff_do);

  /************************************************************/
  /****	THOUGHT OF DOING WITH A SMALLER SHIFT AT STEP 0  ****/
  /****	STILL RECOMBINING WITHIN THE SAME TREE   ,  BUT	 ****/
  /****	ACTUALLY   , THIS METHOD DOES NOT SOLVE THE PB...  ****/
  /************************************************************/

  /****  in order to get a shift in the asset of 0.01	 ****/
  /*************************************************************

          dx_shift_up = log(1.0+0.01/contract.spot)/model.delta_lns;
          dx_shift_do = log(1.0-0.01/contract.spot)/model.delta_lns;

          m = steps[0].drift_h;
          var = steps[0].vol_h;
          m -= steps[0].d_index * model.delta_lns;

          p1_do = (m*m + var - model.delta_lns*m + dx_shift_up*dx_shift_up*
                                  model.delta_lns*model.delta_lns -
                                  dx_shift_up*model.delta_lns*model.delta_lns +
                                  2*m*dx_shift_up*model.delta_lns) /
                                  ( 2*model.delta_lns*model.delta_lns );

          p1_up =(m*m + var + model.delta_lns*m + dx_shift_up*dx_shift_up*
                                  model.delta_lns*model.delta_lns +
                                  dx_shift_up*model.delta_lns*model.delta_lns +
                                  2*m*dx_shift_up*model.delta_lns) /
                                  ( 2*model.delta_lns*model.delta_lns );

          p1_mi = 1 - p1_up - p1_do;

          p2_do = (m*m + var - model.delta_lns*m + dx_shift_do*dx_shift_do*
                                  model.delta_lns*model.delta_lns -
                                  dx_shift_do*model.delta_lns*model.delta_lns +
                                  2*m*dx_shift_do*model.delta_lns) /
                                  ( 2*model.delta_lns*model.delta_lns );

          p2_up =(m*m + var + model.delta_lns*m + dx_shift_do*dx_shift_do*
                                  model.delta_lns*model.delta_lns +
                                  dx_shift_do*model.delta_lns*model.delta_lns +
                                  2*m*dx_shift_do*model.delta_lns) /
                                  ( 2*model.delta_lns*model.delta_lns );

          p2_mi = 1 - p2_up - p2_do;


          price1 = p1_up * prem_prev[1-d0];
          price1 += p1_mi * prem_prev[2-d0];
          price1 += p1_do * prem_prev[3-d0];
          price1 *= disc_fact;

          price2 = p2_up * prem_prev[1-d0];
          price2 += p2_mi * prem_prev[2-d0];
          price2 += p2_do * prem_prev[3-d0];
          price2 *= disc_fact;

          new_delta = (price1-price2)/0.02;

          new_gamma = (price1-price)/0.01;
          new_gamma -= (price-price2)/0.01;
          new_gamma /= 0.01;

  *************************************************************/
  /************************************************************/
  /****	ACTUALLY   , GIVES BACK SAME RESULTS ... HOW NICE  ****/
  /************************************************************/

  theta_1d = (prem_prev[2 - d0] - price) / (model.h * 365.0);

  DEBUG_PRINT &&printf("end of trin \n");

  if (contract.quanto == 1) {
    free_dmatrix(model_fx.date_cumvol, 0, 1, 0, model_fx.size_tsv - 1);
    srt_free(steps_fx);
  }
  free_dmatrix(model.date_cumvol, 0, 1, 0, model.size_tsv - 1);
  srt_free(prem_cur);
  srt_free(prem_prev);
  srt_free(steps);

  BUG_DEBUG &&printf("17 \n");

  if (contract.ret_type == 1)
    return (delta);
  else if (contract.ret_type == 2)
    return (gamma);
  else if (contract.ret_type == 5)
    return (theta_1d);
  else
    return (price);
}

/***************************************************************/
/*********** Payoff and end_tree functions *********************/
/***************************************************************/

double endtree_european(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.intrinsic > 0)
    result = node.intrinsic;
  else
    result = 0;
  return (result);
}

double endtree_partial(Deal contract, Model model, Step step, Node node) {
  double result;

  result =
      srt_f_optblksch(node.forward, contract.strike, step.vol,
                      (double)(contract.opt_mat - contract.sub_period) / 365.0,
                      step.df_swap_end, contract.call_put, 0);
  return (result);
}

double payoff_european(Deal contract, Model model, Step step, Node node) {
  double result;
  if (node.disc_premium > 0)
    result = node.disc_premium;
  else
    result = 0;

  return (result);
}

double payoff_american(Deal contract, Model model, Step step, Node node) {
  double result;
  if (node.disc_premium > node.intrinsic)
    result = node.disc_premium;
  else
    result = node.intrinsic;
  return (result);
}

/***********  SHOUT LOCK  	************/

double endtree_shout_lock(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.intrinsic > 0) {
    result = node.intrinsic * step.df_swap_end;
    result += srt_f_optblksch(node.forward, node.asset, step.vol,
                              (double)(contract.opt_mat - contract.sub_period) /
                                  365.0,
                              step.df_swap_end, contract.call_put, 0);
    /****	replacing step.time_end	by mat - sub 	****/
    /****   because had pb when time_end=0  , no more	****/
  } else
    result = srt_f_optblksch(node.forward, contract.strike, step.vol,
                             (double)(contract.opt_mat - contract.sub_period) /
                                 365.0,
                             step.df_swap_end, contract.call_put, 0);
  return (result);
}

double payoff_shout_lock(Deal contract, Model model, Step step, Node node) {
  double result, prem_val;

  if (node.intrinsic > 0) {
    prem_val = node.intrinsic * step.df_swap_end;
    prem_val +=
        srt_f_optblksch(node.forward, node.asset, step.vol, step.time_end,
                        step.df_swap_end, contract.call_put, 0);
  } else
    prem_val = 0;

  if (node.disc_premium > prem_val)
    result = node.disc_premium;
  else
    result = prem_val;

  return (result);
}

/*********** SHOUT RESET *************/

double endtree_shout_reset(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.intrinsic > 0)
    result = srt_f_optblksch(node.forward, contract.strike, step.vol,
                             (double)(contract.opt_mat - contract.sub_period) /
                                 365.0,
                             step.df_swap_end, contract.call_put, 0);
  else /** RESET **/
    result = srt_f_optblksch(node.forward, node.asset, step.vol,
                             (double)(contract.opt_mat - contract.sub_period) /
                                 365.0,
                             step.df_swap_end, contract.call_put, 0);

  return (result);
}

double payoff_shout_reset(Deal contract, Model model, Step step, Node node) {
  double result, prem_val;

  if (node.intrinsic > 0)
    prem_val = 0;
  else
    prem_val =
        srt_f_optblksch(node.forward, node.asset, step.vol, step.time_end,
                        step.df_swap_end, contract.call_put, 0);

  if (node.disc_premium > prem_val)
    result = node.disc_premium;
  else
    result = prem_val;

  return (result);
}

/*********** SHOUT RESET AND LOCK *************/

double endtree_both_shout(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.intrinsic > 0)

  {
    result = node.intrinsic * step.df_swap_end;
    result += srt_f_optblksch(node.forward, node.asset, step.vol,
                              (double)(contract.opt_mat - contract.sub_period) /
                                  365.0,
                              step.df_swap_end, contract.call_put, 0);
  } else
    result = srt_f_optblksch(node.forward, node.asset, step.vol,
                             (double)(contract.opt_mat - contract.sub_period) /
                                 365.0,
                             step.df_swap_end, contract.call_put, 0);

  return (result);
}

double payoff_both_shout(Deal contract, Model model, Step step, Node node) {
  double result, prem_val;

  if (node.intrinsic > 0) {
    prem_val = node.intrinsic * step.df_swap_end;
    prem_val +=
        srt_f_optblksch(node.forward, node.asset, step.vol, step.time_end,
                        step.df_swap_end, contract.call_put, 0);
  } else
    prem_val =
        srt_f_optblksch(node.forward, node.asset, step.vol, step.time_end,
                        step.df_swap_end, contract.call_put, 0);

  if (node.disc_premium > prem_val)
    result = node.disc_premium;
  else
    result = prem_val;

  return (result);
}

double payoff_instal(Deal contract, Model model, Step step, Node node) {
  double result;

  if (step.flag_date == 1) {
    if (node.disc_premium > step.cpn)
      result = node.disc_premium - step.cpn;
    else
      result = 0;
  } else
    result = node.disc_premium;

  return (result);
}

double endtree_extinguish(Deal contract, Model model, Step step, Node node) {
  double result;

  if ((contract.barrier - node.asset) * (1 - 2 * contract.do_up) >= 0.0)
    result = 0.0;
  else {
    if (contract.sub_period == contract.opt_mat) {
      if (node.intrinsic > 0) {
        if (node.position == model.hit_level) {
          result = ((node.asset + contract.barrier) * 0.5 - contract.strike);
          result *= (contract.barrier - node.asset) /
                    (node.asset * (1 - exp(-model.delta_lns)));
        } else
          result = node.intrinsic;
      } else
        result = 0.0;
    } else
      result = srt_f_optblksch(
          node.forward, contract.strike, step.vol,
          (double)(contract.opt_mat - contract.sub_period) / 365.0,
          step.df_swap_end, contract.call_put, 0);
  }

  return (result);
}

double payoff_extinguish(Deal contract, Model model, Step step, Node node) {
  double result;

  if ((contract.barrier - node.asset) * (1 - 2 * contract.do_up) >= 0.0)
    result = 0.0;
  else {
    if (node.disc_premium > 0)
      result = node.disc_premium;
    else
      result = 0;
  }

  return (result);
}

double payoff_mid_atlantic(Deal contract, Model model, Step step, Node node) {
  double result;

  if (step.flag_date == 1) {
    if (node.disc_premium > node.intrinsic)
      result = node.disc_premium;
    else
      result = node.intrinsic;
  } else
    result = node.disc_premium;

  return (result);
}

/******************** new You Choose Option ***************/

double endtree_you_choose(Deal contract, Model model, Step step, Node node) {
  double result, call_price, put_price;

  call_price =
      srt_f_optblksch(node.forward, contract.strike, step.vol,
                      (double)(contract.opt_mat - contract.sub_period) / 365.0,
                      step.df_swap_end, 0, 0);

  /** strike is call strike **/

  put_price =
      srt_f_optblksch(node.forward, contract.barrier, step.vol,
                      (double)(contract.opt_mat - contract.sub_period) / 365.0,
                      step.df_swap_end, 1, 0);

  /** Barrier is put strike **/

  if (call_price >= put_price) {
    result = call_price;

  } else {
    result = put_price;

    H_DEBUG &&printf(" \n put at =%f ", result);
  }
  return (result);
}

/**********************************************************/
/******************* AMER AND EUROPEAN DIGITALS ***********/
/**********************************************************/

double endtree_digit(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.asset < contract.barrier)

    result = contract.do_up;
  /**** under barrier  , down is out **/

  else
    result = 1 - contract.do_up;
  /**  above barrier  , down pays 1 **/

  /************** not a very good smoothing of payoff cf probas **/

  if ((contract.call_put == 1) && (node.position == model.hit_level)) {
    result = (node.asset - contract.barrier) * (1 - 2 * contract.do_up);
    result +=
        model.diff_do * contract.do_up + model.diff_up * (1 - contract.do_up);
    result /=
        model.diff_do * contract.do_up + model.diff_up * (1 - contract.do_up);

    H_DEBUG &&printf(" \n  new result =%f ", result);
  }

  /*****************************************************************/

  return (result);
}

/***********  digital  paid at the end**********/

double payoff_amer_digit(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.asset < contract.barrier)
    result = node.disc_premium * contract.do_up;

  else
    result = node.disc_premium * (1 - contract.do_up);

  return (result);
}

/*******************************************************/
/*********	  digital immediate pay		********/
/****************** has to be price like a in *********/
/*****************************************************/

double endtree_pay_digit(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.asset > contract.barrier)

    result = contract.do_up;
  /**** above barrier  , up is in **/

  else
    result = 1 - contract.do_up;
  /**  under barrier down is in **/

  return (result);
}

double payoff_pay_digit(Deal contract, Model model, Step step, Node node) {
  double result;

  if (node.asset < contract.barrier) {
    result = node.disc_premium * contract.do_up;
    result += 1 * (1 - contract.do_up);

    /**** do is a down and in ***/
  } else {
    result = 1 * contract.do_up;

    /*** up is a down and in ****/

    result += node.disc_premium * (1 - contract.do_up);
  }

  return (result);
}
