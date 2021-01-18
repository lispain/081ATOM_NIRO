/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_234270915294207647);
void inv_err_fun(double *nom_x, double *true_x, double *out_3520324582001283245);
void H_mod_fun(double *state, double *out_2988122354651275992);
void f_fun(double *state, double dt, double *out_7779891258215356160);
void F_fun(double *state, double dt, double *out_7589799868175156848);
void h_3(double *state, double *unused, double *out_5250311108539441373);
void H_3(double *state, double *unused, double *out_7869960643351505062);
void h_4(double *state, double *unused, double *out_5206639634793134209);
void H_4(double *state, double *unused, double *out_8256083701102277576);
void h_9(double *state, double *unused, double *out_2519245694720930710);
void H_9(double *state, double *unused, double *out_5933578157231421606);
void h_10(double *state, double *unused, double *out_4759890688823990129);
void H_10(double *state, double *unused, double *out_4027739071800323038);
void h_12(double *state, double *unused, double *out_7182999656089221228);
void H_12(double *state, double *unused, double *out_6149528189363068007);
void h_31(double *state, double *unused, double *out_6638058017871132724);
void H_31(double *state, double *unused, double *out_8824617964506487828);
void h_32(double *state, double *unused, double *out_7294741787366483558);
void H_32(double *state, double *unused, double *out_706622715902999258);
void h_13(double *state, double *unused, double *out_8453116867543887942);
void H_13(double *state, double *unused, double *out_724568766321066001);
void h_14(double *state, double *unused, double *out_2519245694720930710);
void H_14(double *state, double *unused, double *out_5933578157231421606);
void h_19(double *state, double *unused, double *out_8019342888595656710);
void H_19(double *state, double *unused, double *out_1500898306128990425);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);