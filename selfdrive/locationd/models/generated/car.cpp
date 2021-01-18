
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4159821151823311033) {
   out_4159821151823311033[0] = delta_x[0] + nom_x[0];
   out_4159821151823311033[1] = delta_x[1] + nom_x[1];
   out_4159821151823311033[2] = delta_x[2] + nom_x[2];
   out_4159821151823311033[3] = delta_x[3] + nom_x[3];
   out_4159821151823311033[4] = delta_x[4] + nom_x[4];
   out_4159821151823311033[5] = delta_x[5] + nom_x[5];
   out_4159821151823311033[6] = delta_x[6] + nom_x[6];
   out_4159821151823311033[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4974507302638080604) {
   out_4974507302638080604[0] = -nom_x[0] + true_x[0];
   out_4974507302638080604[1] = -nom_x[1] + true_x[1];
   out_4974507302638080604[2] = -nom_x[2] + true_x[2];
   out_4974507302638080604[3] = -nom_x[3] + true_x[3];
   out_4974507302638080604[4] = -nom_x[4] + true_x[4];
   out_4974507302638080604[5] = -nom_x[5] + true_x[5];
   out_4974507302638080604[6] = -nom_x[6] + true_x[6];
   out_4974507302638080604[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_6763524490318930976) {
   out_6763524490318930976[0] = 1.0;
   out_6763524490318930976[1] = 0.0;
   out_6763524490318930976[2] = 0.0;
   out_6763524490318930976[3] = 0.0;
   out_6763524490318930976[4] = 0.0;
   out_6763524490318930976[5] = 0.0;
   out_6763524490318930976[6] = 0.0;
   out_6763524490318930976[7] = 0.0;
   out_6763524490318930976[8] = 0.0;
   out_6763524490318930976[9] = 1.0;
   out_6763524490318930976[10] = 0.0;
   out_6763524490318930976[11] = 0.0;
   out_6763524490318930976[12] = 0.0;
   out_6763524490318930976[13] = 0.0;
   out_6763524490318930976[14] = 0.0;
   out_6763524490318930976[15] = 0.0;
   out_6763524490318930976[16] = 0.0;
   out_6763524490318930976[17] = 0.0;
   out_6763524490318930976[18] = 1.0;
   out_6763524490318930976[19] = 0.0;
   out_6763524490318930976[20] = 0.0;
   out_6763524490318930976[21] = 0.0;
   out_6763524490318930976[22] = 0.0;
   out_6763524490318930976[23] = 0.0;
   out_6763524490318930976[24] = 0.0;
   out_6763524490318930976[25] = 0.0;
   out_6763524490318930976[26] = 0.0;
   out_6763524490318930976[27] = 1.0;
   out_6763524490318930976[28] = 0.0;
   out_6763524490318930976[29] = 0.0;
   out_6763524490318930976[30] = 0.0;
   out_6763524490318930976[31] = 0.0;
   out_6763524490318930976[32] = 0.0;
   out_6763524490318930976[33] = 0.0;
   out_6763524490318930976[34] = 0.0;
   out_6763524490318930976[35] = 0.0;
   out_6763524490318930976[36] = 1.0;
   out_6763524490318930976[37] = 0.0;
   out_6763524490318930976[38] = 0.0;
   out_6763524490318930976[39] = 0.0;
   out_6763524490318930976[40] = 0.0;
   out_6763524490318930976[41] = 0.0;
   out_6763524490318930976[42] = 0.0;
   out_6763524490318930976[43] = 0.0;
   out_6763524490318930976[44] = 0.0;
   out_6763524490318930976[45] = 1.0;
   out_6763524490318930976[46] = 0.0;
   out_6763524490318930976[47] = 0.0;
   out_6763524490318930976[48] = 0.0;
   out_6763524490318930976[49] = 0.0;
   out_6763524490318930976[50] = 0.0;
   out_6763524490318930976[51] = 0.0;
   out_6763524490318930976[52] = 0.0;
   out_6763524490318930976[53] = 0.0;
   out_6763524490318930976[54] = 1.0;
   out_6763524490318930976[55] = 0.0;
   out_6763524490318930976[56] = 0.0;
   out_6763524490318930976[57] = 0.0;
   out_6763524490318930976[58] = 0.0;
   out_6763524490318930976[59] = 0.0;
   out_6763524490318930976[60] = 0.0;
   out_6763524490318930976[61] = 0.0;
   out_6763524490318930976[62] = 0.0;
   out_6763524490318930976[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9114670843477850011) {
   out_9114670843477850011[0] = state[0];
   out_9114670843477850011[1] = state[1];
   out_9114670843477850011[2] = state[2];
   out_9114670843477850011[3] = state[3];
   out_9114670843477850011[4] = state[4];
   out_9114670843477850011[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9114670843477850011[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9114670843477850011[7] = state[7];
}
void F_fun(double *state, double dt, double *out_6681288621460938018) {
   out_6681288621460938018[0] = 1;
   out_6681288621460938018[1] = 0;
   out_6681288621460938018[2] = 0;
   out_6681288621460938018[3] = 0;
   out_6681288621460938018[4] = 0;
   out_6681288621460938018[5] = 0;
   out_6681288621460938018[6] = 0;
   out_6681288621460938018[7] = 0;
   out_6681288621460938018[8] = 0;
   out_6681288621460938018[9] = 1;
   out_6681288621460938018[10] = 0;
   out_6681288621460938018[11] = 0;
   out_6681288621460938018[12] = 0;
   out_6681288621460938018[13] = 0;
   out_6681288621460938018[14] = 0;
   out_6681288621460938018[15] = 0;
   out_6681288621460938018[16] = 0;
   out_6681288621460938018[17] = 0;
   out_6681288621460938018[18] = 1;
   out_6681288621460938018[19] = 0;
   out_6681288621460938018[20] = 0;
   out_6681288621460938018[21] = 0;
   out_6681288621460938018[22] = 0;
   out_6681288621460938018[23] = 0;
   out_6681288621460938018[24] = 0;
   out_6681288621460938018[25] = 0;
   out_6681288621460938018[26] = 0;
   out_6681288621460938018[27] = 1;
   out_6681288621460938018[28] = 0;
   out_6681288621460938018[29] = 0;
   out_6681288621460938018[30] = 0;
   out_6681288621460938018[31] = 0;
   out_6681288621460938018[32] = 0;
   out_6681288621460938018[33] = 0;
   out_6681288621460938018[34] = 0;
   out_6681288621460938018[35] = 0;
   out_6681288621460938018[36] = 1;
   out_6681288621460938018[37] = 0;
   out_6681288621460938018[38] = 0;
   out_6681288621460938018[39] = 0;
   out_6681288621460938018[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6681288621460938018[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6681288621460938018[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6681288621460938018[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6681288621460938018[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6681288621460938018[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6681288621460938018[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6681288621460938018[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6681288621460938018[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6681288621460938018[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6681288621460938018[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6681288621460938018[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6681288621460938018[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6681288621460938018[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6681288621460938018[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6681288621460938018[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6681288621460938018[56] = 0;
   out_6681288621460938018[57] = 0;
   out_6681288621460938018[58] = 0;
   out_6681288621460938018[59] = 0;
   out_6681288621460938018[60] = 0;
   out_6681288621460938018[61] = 0;
   out_6681288621460938018[62] = 0;
   out_6681288621460938018[63] = 1;
}
void h_25(double *state, double *unused, double *out_104949727670378771) {
   out_104949727670378771[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1385266769028112500) {
   out_1385266769028112500[0] = 0;
   out_1385266769028112500[1] = 0;
   out_1385266769028112500[2] = 0;
   out_1385266769028112500[3] = 0;
   out_1385266769028112500[4] = 0;
   out_1385266769028112500[5] = 0;
   out_1385266769028112500[6] = 1;
   out_1385266769028112500[7] = 0;
}
void h_24(double *state, double *unused, double *out_6306333318975887828) {
   out_6306333318975887828[0] = state[4];
   out_6306333318975887828[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3082919364583970070) {
   out_3082919364583970070[0] = 0;
   out_3082919364583970070[1] = 0;
   out_3082919364583970070[2] = 0;
   out_3082919364583970070[3] = 0;
   out_3082919364583970070[4] = 1;
   out_3082919364583970070[5] = 0;
   out_3082919364583970070[6] = 0;
   out_3082919364583970070[7] = 0;
   out_3082919364583970070[8] = 0;
   out_3082919364583970070[9] = 0;
   out_3082919364583970070[10] = 0;
   out_3082919364583970070[11] = 0;
   out_3082919364583970070[12] = 0;
   out_3082919364583970070[13] = 1;
   out_3082919364583970070[14] = 0;
   out_3082919364583970070[15] = 0;
}
void h_30(double *state, double *unused, double *out_8102547018412528079) {
   out_8102547018412528079[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8228632141921070780) {
   out_8228632141921070780[0] = 0;
   out_8228632141921070780[1] = 0;
   out_8228632141921070780[2] = 0;
   out_8228632141921070780[3] = 0;
   out_8228632141921070780[4] = 1;
   out_8228632141921070780[5] = 0;
   out_8228632141921070780[6] = 0;
   out_8228632141921070780[7] = 0;
}
void h_26(double *state, double *unused, double *out_6412099183676847119) {
   out_6412099183676847119[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2691502540982953645) {
   out_2691502540982953645[0] = 0;
   out_2691502540982953645[1] = 0;
   out_2691502540982953645[2] = 0;
   out_2691502540982953645[3] = 0;
   out_2691502540982953645[4] = 0;
   out_2691502540982953645[5] = 0;
   out_2691502540982953645[6] = 0;
   out_2691502540982953645[7] = 1;
}
void h_27(double *state, double *unused, double *out_2499455593138318880) {
   out_2499455593138318880[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8930529943951855524) {
   out_8930529943951855524[0] = 0;
   out_8930529943951855524[1] = 0;
   out_8930529943951855524[2] = 0;
   out_8930529943951855524[3] = 1;
   out_8930529943951855524[4] = 0;
   out_8930529943951855524[5] = 0;
   out_8930529943951855524[6] = 0;
   out_8930529943951855524[7] = 0;
}
void h_29(double *state, double *unused, double *out_1362227883290182772) {
   out_1362227883290182772[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6394313162659648883) {
   out_6394313162659648883[0] = 0;
   out_6394313162659648883[1] = 1;
   out_6394313162659648883[2] = 0;
   out_6394313162659648883[3] = 0;
   out_6394313162659648883[4] = 0;
   out_6394313162659648883[5] = 0;
   out_6394313162659648883[6] = 0;
   out_6394313162659648883[7] = 0;
}
void h_28(double *state, double *unused, double *out_7326742220570034943) {
   out_7326742220570034943[0] = state[5];
   out_7326742220570034943[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3934913113652524796) {
   out_3934913113652524796[0] = 0;
   out_3934913113652524796[1] = 0;
   out_3934913113652524796[2] = 0;
   out_3934913113652524796[3] = 0;
   out_3934913113652524796[4] = 0;
   out_3934913113652524796[5] = 1;
   out_3934913113652524796[6] = 0;
   out_3934913113652524796[7] = 0;
   out_3934913113652524796[8] = 0;
   out_3934913113652524796[9] = 0;
   out_3934913113652524796[10] = 0;
   out_3934913113652524796[11] = 0;
   out_3934913113652524796[12] = 0;
   out_3934913113652524796[13] = 0;
   out_3934913113652524796[14] = 1;
   out_3934913113652524796[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
