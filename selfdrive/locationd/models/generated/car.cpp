
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
void err_fun(double *nom_x, double *delta_x, double *out_718929463194776553) {
   out_718929463194776553[0] = delta_x[0] + nom_x[0];
   out_718929463194776553[1] = delta_x[1] + nom_x[1];
   out_718929463194776553[2] = delta_x[2] + nom_x[2];
   out_718929463194776553[3] = delta_x[3] + nom_x[3];
   out_718929463194776553[4] = delta_x[4] + nom_x[4];
   out_718929463194776553[5] = delta_x[5] + nom_x[5];
   out_718929463194776553[6] = delta_x[6] + nom_x[6];
   out_718929463194776553[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3692640817737576391) {
   out_3692640817737576391[0] = -nom_x[0] + true_x[0];
   out_3692640817737576391[1] = -nom_x[1] + true_x[1];
   out_3692640817737576391[2] = -nom_x[2] + true_x[2];
   out_3692640817737576391[3] = -nom_x[3] + true_x[3];
   out_3692640817737576391[4] = -nom_x[4] + true_x[4];
   out_3692640817737576391[5] = -nom_x[5] + true_x[5];
   out_3692640817737576391[6] = -nom_x[6] + true_x[6];
   out_3692640817737576391[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_8200148637365978873) {
   out_8200148637365978873[0] = 1.0;
   out_8200148637365978873[1] = 0.0;
   out_8200148637365978873[2] = 0.0;
   out_8200148637365978873[3] = 0.0;
   out_8200148637365978873[4] = 0.0;
   out_8200148637365978873[5] = 0.0;
   out_8200148637365978873[6] = 0.0;
   out_8200148637365978873[7] = 0.0;
   out_8200148637365978873[8] = 0.0;
   out_8200148637365978873[9] = 1.0;
   out_8200148637365978873[10] = 0.0;
   out_8200148637365978873[11] = 0.0;
   out_8200148637365978873[12] = 0.0;
   out_8200148637365978873[13] = 0.0;
   out_8200148637365978873[14] = 0.0;
   out_8200148637365978873[15] = 0.0;
   out_8200148637365978873[16] = 0.0;
   out_8200148637365978873[17] = 0.0;
   out_8200148637365978873[18] = 1.0;
   out_8200148637365978873[19] = 0.0;
   out_8200148637365978873[20] = 0.0;
   out_8200148637365978873[21] = 0.0;
   out_8200148637365978873[22] = 0.0;
   out_8200148637365978873[23] = 0.0;
   out_8200148637365978873[24] = 0.0;
   out_8200148637365978873[25] = 0.0;
   out_8200148637365978873[26] = 0.0;
   out_8200148637365978873[27] = 1.0;
   out_8200148637365978873[28] = 0.0;
   out_8200148637365978873[29] = 0.0;
   out_8200148637365978873[30] = 0.0;
   out_8200148637365978873[31] = 0.0;
   out_8200148637365978873[32] = 0.0;
   out_8200148637365978873[33] = 0.0;
   out_8200148637365978873[34] = 0.0;
   out_8200148637365978873[35] = 0.0;
   out_8200148637365978873[36] = 1.0;
   out_8200148637365978873[37] = 0.0;
   out_8200148637365978873[38] = 0.0;
   out_8200148637365978873[39] = 0.0;
   out_8200148637365978873[40] = 0.0;
   out_8200148637365978873[41] = 0.0;
   out_8200148637365978873[42] = 0.0;
   out_8200148637365978873[43] = 0.0;
   out_8200148637365978873[44] = 0.0;
   out_8200148637365978873[45] = 1.0;
   out_8200148637365978873[46] = 0.0;
   out_8200148637365978873[47] = 0.0;
   out_8200148637365978873[48] = 0.0;
   out_8200148637365978873[49] = 0.0;
   out_8200148637365978873[50] = 0.0;
   out_8200148637365978873[51] = 0.0;
   out_8200148637365978873[52] = 0.0;
   out_8200148637365978873[53] = 0.0;
   out_8200148637365978873[54] = 1.0;
   out_8200148637365978873[55] = 0.0;
   out_8200148637365978873[56] = 0.0;
   out_8200148637365978873[57] = 0.0;
   out_8200148637365978873[58] = 0.0;
   out_8200148637365978873[59] = 0.0;
   out_8200148637365978873[60] = 0.0;
   out_8200148637365978873[61] = 0.0;
   out_8200148637365978873[62] = 0.0;
   out_8200148637365978873[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8285677604977850576) {
   out_8285677604977850576[0] = state[0];
   out_8285677604977850576[1] = state[1];
   out_8285677604977850576[2] = state[2];
   out_8285677604977850576[3] = state[3];
   out_8285677604977850576[4] = state[4];
   out_8285677604977850576[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8285677604977850576[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8285677604977850576[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8095088944597024994) {
   out_8095088944597024994[0] = 1;
   out_8095088944597024994[1] = 0;
   out_8095088944597024994[2] = 0;
   out_8095088944597024994[3] = 0;
   out_8095088944597024994[4] = 0;
   out_8095088944597024994[5] = 0;
   out_8095088944597024994[6] = 0;
   out_8095088944597024994[7] = 0;
   out_8095088944597024994[8] = 0;
   out_8095088944597024994[9] = 1;
   out_8095088944597024994[10] = 0;
   out_8095088944597024994[11] = 0;
   out_8095088944597024994[12] = 0;
   out_8095088944597024994[13] = 0;
   out_8095088944597024994[14] = 0;
   out_8095088944597024994[15] = 0;
   out_8095088944597024994[16] = 0;
   out_8095088944597024994[17] = 0;
   out_8095088944597024994[18] = 1;
   out_8095088944597024994[19] = 0;
   out_8095088944597024994[20] = 0;
   out_8095088944597024994[21] = 0;
   out_8095088944597024994[22] = 0;
   out_8095088944597024994[23] = 0;
   out_8095088944597024994[24] = 0;
   out_8095088944597024994[25] = 0;
   out_8095088944597024994[26] = 0;
   out_8095088944597024994[27] = 1;
   out_8095088944597024994[28] = 0;
   out_8095088944597024994[29] = 0;
   out_8095088944597024994[30] = 0;
   out_8095088944597024994[31] = 0;
   out_8095088944597024994[32] = 0;
   out_8095088944597024994[33] = 0;
   out_8095088944597024994[34] = 0;
   out_8095088944597024994[35] = 0;
   out_8095088944597024994[36] = 1;
   out_8095088944597024994[37] = 0;
   out_8095088944597024994[38] = 0;
   out_8095088944597024994[39] = 0;
   out_8095088944597024994[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8095088944597024994[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8095088944597024994[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8095088944597024994[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8095088944597024994[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8095088944597024994[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8095088944597024994[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8095088944597024994[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8095088944597024994[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8095088944597024994[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8095088944597024994[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8095088944597024994[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8095088944597024994[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8095088944597024994[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8095088944597024994[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8095088944597024994[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8095088944597024994[56] = 0;
   out_8095088944597024994[57] = 0;
   out_8095088944597024994[58] = 0;
   out_8095088944597024994[59] = 0;
   out_8095088944597024994[60] = 0;
   out_8095088944597024994[61] = 0;
   out_8095088944597024994[62] = 0;
   out_8095088944597024994[63] = 1;
}
void h_25(double *state, double *unused, double *out_7766190187124533728) {
   out_7766190187124533728[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6634415785610795242) {
   out_6634415785610795242[0] = 0;
   out_6634415785610795242[1] = 0;
   out_6634415785610795242[2] = 0;
   out_6634415785610795242[3] = 0;
   out_6634415785610795242[4] = 0;
   out_6634415785610795242[5] = 0;
   out_6634415785610795242[6] = 1;
   out_6634415785610795242[7] = 0;
}
void h_24(double *state, double *unused, double *out_1129534900983040369) {
   out_1129534900983040369[0] = state[4];
   out_1129534900983040369[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7896804849806535374) {
   out_7896804849806535374[0] = 0;
   out_7896804849806535374[1] = 0;
   out_7896804849806535374[2] = 0;
   out_7896804849806535374[3] = 0;
   out_7896804849806535374[4] = 1;
   out_7896804849806535374[5] = 0;
   out_7896804849806535374[6] = 0;
   out_7896804849806535374[7] = 0;
   out_7896804849806535374[8] = 0;
   out_7896804849806535374[9] = 0;
   out_7896804849806535374[10] = 0;
   out_7896804849806535374[11] = 0;
   out_7896804849806535374[12] = 0;
   out_7896804849806535374[13] = 1;
   out_7896804849806535374[14] = 0;
   out_7896804849806535374[15] = 0;
}
void h_30(double *state, double *unused, double *out_115082463629515457) {
   out_115082463629515457[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2979483125338388038) {
   out_2979483125338388038[0] = 0;
   out_2979483125338388038[1] = 0;
   out_2979483125338388038[2] = 0;
   out_2979483125338388038[3] = 0;
   out_2979483125338388038[4] = 1;
   out_2979483125338388038[5] = 0;
   out_2979483125338388038[6] = 0;
   out_2979483125338388038[7] = 0;
}
void h_26(double *state, double *unused, double *out_488233421046962060) {
   out_488233421046962060[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5292979651915147690) {
   out_5292979651915147690[0] = 0;
   out_5292979651915147690[1] = 0;
   out_5292979651915147690[2] = 0;
   out_5292979651915147690[3] = 0;
   out_5292979651915147690[4] = 0;
   out_5292979651915147690[5] = 0;
   out_5292979651915147690[6] = 0;
   out_5292979651915147690[7] = 1;
}
void h_27(double *state, double *unused, double *out_8709700040026196163) {
   out_8709700040026196163[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4267065113175013350) {
   out_4267065113175013350[0] = 0;
   out_4267065113175013350[1] = 0;
   out_4267065113175013350[2] = 0;
   out_4267065113175013350[3] = 1;
   out_4267065113175013350[4] = 0;
   out_4267065113175013350[5] = 0;
   out_4267065113175013350[6] = 0;
   out_4267065113175013350[7] = 0;
}
void h_29(double *state, double *unused, double *out_828427389272146978) {
   out_828427389272146978[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6803281894467219991) {
   out_6803281894467219991[0] = 0;
   out_6803281894467219991[1] = 1;
   out_6803281894467219991[2] = 0;
   out_6803281894467219991[3] = 0;
   out_6803281894467219991[4] = 0;
   out_6803281894467219991[5] = 0;
   out_6803281894467219991[6] = 0;
   out_6803281894467219991[7] = 0;
}
void h_28(double *state, double *unused, double *out_4420494587383410260) {
   out_4420494587383410260[0] = state[5];
   out_4420494587383410260[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3532106745666521376) {
   out_3532106745666521376[0] = 0;
   out_3532106745666521376[1] = 0;
   out_3532106745666521376[2] = 0;
   out_3532106745666521376[3] = 0;
   out_3532106745666521376[4] = 0;
   out_3532106745666521376[5] = 1;
   out_3532106745666521376[6] = 0;
   out_3532106745666521376[7] = 0;
   out_3532106745666521376[8] = 0;
   out_3532106745666521376[9] = 0;
   out_3532106745666521376[10] = 0;
   out_3532106745666521376[11] = 0;
   out_3532106745666521376[12] = 0;
   out_3532106745666521376[13] = 0;
   out_3532106745666521376[14] = 1;
   out_3532106745666521376[15] = 0;
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
