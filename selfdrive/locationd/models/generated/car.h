/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4159821151823311033);
void inv_err_fun(double *nom_x, double *true_x, double *out_4974507302638080604);
void H_mod_fun(double *state, double *out_6763524490318930976);
void f_fun(double *state, double dt, double *out_9114670843477850011);
void F_fun(double *state, double dt, double *out_6681288621460938018);
void h_25(double *state, double *unused, double *out_104949727670378771);
void H_25(double *state, double *unused, double *out_1385266769028112500);
void h_24(double *state, double *unused, double *out_6306333318975887828);
void H_24(double *state, double *unused, double *out_3082919364583970070);
void h_30(double *state, double *unused, double *out_8102547018412528079);
void H_30(double *state, double *unused, double *out_8228632141921070780);
void h_26(double *state, double *unused, double *out_6412099183676847119);
void H_26(double *state, double *unused, double *out_2691502540982953645);
void h_27(double *state, double *unused, double *out_2499455593138318880);
void H_27(double *state, double *unused, double *out_8930529943951855524);
void h_29(double *state, double *unused, double *out_1362227883290182772);
void H_29(double *state, double *unused, double *out_6394313162659648883);
void h_28(double *state, double *unused, double *out_7326742220570034943);
void H_28(double *state, double *unused, double *out_3934913113652524796);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
