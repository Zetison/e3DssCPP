// $Id$
//==============================================================================
//!
//! \file e3Dss.C
//!
//! \date Sep 13 2019
//!
//! \author Jon Vegard Venås / SINTEF
//!
//! \brief Class for exact 3D scattering solutions.
//!
//==============================================================================



#ifndef _E3DSS_H
#define _E3DSS_H

#include <iostream>
#include <complex>
#include <vector>

typedef double Real;
enum neumannConditions {PLANE_WAVE, POINT_SOURCE_WAVE};
enum internalBC {SHBC, SSBC, ESBC, NNBC};
struct sc
{
  Real r;
  Real theta;
};

class E3Dss
{
public:
  //! \brief The default constructor initializes a rigid spherical shell problem
  E3Dss();
  void set_r_s(Real);
private:
  // Physical properties (constant)
  Vec3 d_vec[3]; //!< Direction of the incident wave
  Real omega; //!< Angular frequency
  std::vector<Real> R_i; //!< Inner radii
  std::vector<Real> R_o; //!< Outer radii
  std::vector<Real> E; //!< Youngs modulus for solid layers
  std::vector<Real> nu; //!< Poisson ratio for solid layers
  std::vector<Real> rho_s; //!< Mass densities of solid layers
  std::vector<Real> rho_f; //!< Mass densities of fluid layers
  std::vector<Real> c_f; //!< Speed of sound in fluid layers
  Real Eps; //!< Small parameter for series truncation
  size_t N_max; //!< Upper limit for the number of terms in the series
  bool calc_farField; //!< Set true if the far field is to be calculated for outer fluid domain
  int neumannCond; //!< 1 = PLANE_WAVE, 2 = POINT_SOURCE_WAVE
  int iBC; //! 1 = SHBC, 2 = SSBC, 3 = ESBC, 4 = NNBC
  Real r_s;	//!< Radius to source location for point charge incident waves
  bool calc_p; //!< Calculate the scattered pressure
  bool calc_dpdx; //!< Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate x
  bool calc_dpdy; //!< Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate y
  bool calc_dpdz; //!< Calculate the derivative of the scattered pressure fields w.r.t. the cartesian coordinate z
  bool calc_u_x; //!< Calculate the x component of the displacement
  bool calc_u_y; //!< Calculate the y component of the displacement
  bool calc_u_z; //!< Calculate the z component of the displacement
  bool calc_du_xdx; //!< Calculate the derivative of the x component of the displacement field w.r.t. x
  bool calc_du_xdy; //!< Calculate the derivative of the x component of the displacement field w.r.t. y
  bool calc_du_xdz; //!< Calculate the derivative of the x component of the displacement field w.r.t. z
  bool calc_du_ydx; //!< Calculate the derivative of the y component of the displacement field w.r.t. x
 bool  calc_du_ydy; //!< Calculate the derivative of the y component of the displacement field w.r.t. y
  bool calc_du_ydz; //!< Calculate the derivative of the y component of the displacement field w.r.t. z
  bool calc_du_zdx; //!< Calculate the derivative of the z component of the displacement field w.r.t. x
 bool  calc_du_zdy; //!< Calculate the derivative of the z component of the displacement field w.r.t. y
  bool calc_du_zdz; //!< Calculate the derivative of the z component of the displacement field w.r.t. z
  bool calc_sigma_xx; //!< Calculate the derivative of the xx component of the stress field (cartesian coordinates)
  bool calc_sigma_yy; //!< Calculate the derivative of the yy component of the stress field (cartesian coordinates)
  bool calc_sigma_zz; //!< Calculate the derivative of the zz component of the stress field (cartesian coordinates)
  bool calc_sigma_yz; //!< Calculate the derivative of the yz component of the stress field (cartesian coordinates)
  bool calc_sigma_xz; //!< Calculate the derivative of the xz component of the stress field (cartesian coordinates)
  bool calc_sigma_xy; //!< Calculate the derivative of the xy component of the stress field (cartesian coordinates)
  bool calc_p_laplace; //!< Calculate the Laplace operator of the scattered pressure fields
  bool calc_errorsDisplacementCondition; //!< Calculate the errors for the displacement conditions
  bool calc_errorsPressureCondition; //!< Calculate the errors for the pressure conditions
  bool calc_errorsHelmholtz; //!< Calculate errors for the Helmholtz equation
  bool calc_errorsNavier;	//!< Calculate errors for the Navier equation
};

