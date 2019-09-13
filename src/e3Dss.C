// $Id$
//==============================================================================
//!
//! \file e3Dss.h
//!
//! \date Sep 13 2019
//!
//! \author Jon Vegard Venås / SINTEF
//!
//! \brief Class for exact 3D scattering solutions.
//!
//==============================================================================

#include <math.h>
#include "e3Dss.h"
#include "mkl.h"

E3Dss::E3Dss()
{
  d_vec = Vec3({0.0,0.0,1.0});
  omega = 2*pi*1e3
  R_o.push_back(1);
  c_f.push_back(1500);
  Eps = DBL_EPSILON;
  N_max = INFINITY;
  calc_farField = false;
  neumannCond = 1;
  calc_p = true;  
  calc_dpdx = false;
  calc_dpdy = false;
  calc_dpdz = false;
  calc_u_x = false;
  calc_u_y = false;
  calc_u_z = false;
  calc_du_xdx = false;
  calc_du_xdy = false;
  calc_du_xdz = false;
  calc_du_ydx = false;
  calc_du_ydy = false;
  calc_du_ydz = false;
  calc_du_zdx = false;
  calc_du_zdy = false;
  calc_du_zdz = false;
  calc_sigma_xx = false;
  calc_sigma_yy = false;
  calc_sigma_zz = false;
  calc_sigma_yz = false;
  calc_sigma_xz = false;
  calc_sigma_xy = false;
  calc_p_laplace = false;
  calc_errorsDisplacementCondition = false;
  calc_errorsPressureCondition = false;
  calc_errorsHelmholtz= false;    
  calc_errorsNavier = false;
  if(E.size() < R_o.size())
    iBC = SHBC;
  else if(R_o.size() > R_i.size())
    iBC = ESBC;
  else if(R_o.size() == rho_f.size())
    iBC = SSBC;
  else
    iBC = NNBC;
  d_vec.normalize(eps);
}
void E3Dss::set_r_s(Real r_s)
{
  if(r_s <= this->R_0[0]){
    fprintf(stderr, "r_s must satisfy |r_s| > R_o[0]");
    exit(-1);
  }
  this->r_s = r_s;
  std::cout << "'It is not implemented an efficient routine to evaluate Equation (D.7).\n";
}
 E3Dss::computeFields(real **X)
{
  Vec3 e1(0,0,0);

