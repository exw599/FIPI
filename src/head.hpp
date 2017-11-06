/* ----------------------------------------------------------------------
 *
 *   	FIPI - Fast Interface Particle Interactions
 *
 *	https://bottogroup.wordpress.com/
 *
 *   	Chuan Gu, c.gu@qmul.ac.uk
 *
 *   	Copyright (2016) Botto Research Group
 *
 *   	This software is distributed under the GNU General Public License.
 *
 *------------------------------------------------------------------------- */

#ifndef HEAD 
#define HEAD 

#include<fftw3.h>
#include<cmath>
#include<ctime>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<cstring>

const double pi = 3.1415926;

const int Nx = 64;  
const int Ny = 64;  
const int Nz = 64;  

const double Lx = 2*pi;
const double Ly = 2*pi;
const double Lz = 2*pi;

const double dx = Lx/Nx;
const double dy = Ly/Ny;
const double dz = Lz/Nz;

const double dt = 1e-3; // time step
const int t_end = 2e4;  // number of time step

const double Bo = 0.0;
const double sigma   = 0.07;
const double epsilon = 3.0*dx/pow(2,0.5)/2.0; // half width = 3*dx 
const double lambda  = sigma*epsilon*3.0/2.0/pow(2,0.5); // gradient energy coefficient  

const double r_drop = 0.4*Lx;

const double gravity = -1e1;// gravitational acceleration
const double rho1  = 1e-3; //density of fluid
const double rho2  = rho1 - Bo*sigma/gravity/r_drop/r_drop;
const double mu = 1e-3;   //dynamic viscosity of fluid
 
//----------------------------------------------------------------
const double A = 2; // constant in front of attachment force F_pi
const double mobility = pow(0.2*epsilon,2)/mu;
const double r_p = 0.03*r_drop; //solid particle radius
const double rho_p = rho1;
//----------------------------------------------------------------

const int N_p = 20000;
const double r_eq = 5*r_p;
const double r_c  = 5*r_p;
const double r_t  = r_p;

const double kt = 0.0; // thermal energy

#include"parameter.hpp"

const bool switch_particle = true;
const bool switch_ch = true;
const bool switch_ns = true;
#endif
