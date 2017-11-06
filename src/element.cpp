/*----------------------------------------------------------------------
 *  
 *  	FIPI - Fast Interface Particle Interactions
 *   
 *     	https://bottogroup.wordpress.com/
 *     
 *     	Chuan Gu, c.gu@qmul.ac.uk
 *    
 *    	Copyright (2016) Botto Research Group
 *   
 * 	This software is distributed under the GNU General Public License.
 *  
 * ------------------------------------------------------------------------- */

#include"head.hpp"
#include"primitive.hpp"
#include"element.hpp"

//-----------------------------------------------------------------CLASS ELEMENT-------------------------------------------------------	

//****************
element::element()
//****************
{
        r   = r_p;
        
	rho = rho_p;

	flag = 0;

	u     = vector(0.0,0.0,0.0);
	u_old = vector(0.0,0.0,0.0);
        
	c = 1.0;
}

//********************************************
void element::get_local_u(const field& u_grid)
//********************************************
{
	uf.set_value(0.0,0.0,0.0);

	//! use interpolation with wider support
	/*	
	uf.set_value
	(
	u_grid.interpolate(p, 0),
	u_grid.interpolate(p, 1),
	u_grid.interpolate(p, 2)
	);
	*/
}

//********************************************
void element::get_local_c(const field& c_grid)
//********************************************
{
 	c = c_grid.interpolate(p, 0);
}

//*********************************************
void element::get_local_dcdx(const field& grad_c_grid)
//*********************************************
{
	dcdx.set_value
	(	
	grad_c_grid.interpolate(p, 0),
	grad_c_grid.interpolate(p, 1),
	grad_c_grid.interpolate(p, 2)
	);
}


//***********************
void element::get_f_pi()
//***********************
{
        double f_mag = 0.0; // magnitude of adhesion force between interface and a particle
        double d = 0.0;     // distance of particle center to interface

	//!
	double phi = c;  	       
 
        if ( std::abs(phi) <= tanh(2.0) )
        {
	if (dcdx.norm() != 0.0)
	{
		d = 0.5*log((1.0+phi)/(1.0-phi))*sqrt(2.0)*epsilon;	
		if (std::abs(d)<r) f_mag = -A*pi*d*sigma; 
		else	           f_mag = 0.0;
		f_pi = dcdx/dcdx.norm()*f_mag;					
	}
	else	f_pi.set_value(0.0,0.0,0.0); 
        }
        else	f_pi.set_value(0.0,0.0,0.0); 
	
}

//****************************
void element::get_f_ext(int t)
//****************************
{
	f_ext.set_value(0.0,0.0,0.0); 	
}



//**************************************
void element::get_f_rand()
//**************************************
{
        double r1, r2;
        double g1, g2, g3;

        g1 = 11.0;
        g2 = 11.0;
        g3 = 11.0;

        do {
                r1 = double(rand() % 1000 )/1000;
                r2 = double(rand() % 1000 )/1000;
                g1 = sqrt(-2*log(r1))*cos(2*pi*r2);
        } while (std::abs(g1) > 5.0);

        do {
                r1 = double(rand() % 1000 )/1000;
                r2 = double(rand() % 1000 )/1000;
                g2 = sqrt(-2*log(r1))*cos(2*pi*r2);
        } while (std::abs(g2) > 5.0);

        do {
                r1 = double(rand() % 1000 )/1000;
                r2 = double(rand() % 1000 )/1000;
                g3 = sqrt(-2*log(r1))*cos(2*pi*r2);
        } while (std::abs(g3) > 5.0);


        f_rand.set_value
	(
		sqrt(12.0*pi*r_p*mu*kt/dt)*g1,
		sqrt(12.0*pi*r_p*mu*kt/dt)*g2,
		sqrt(12.0*pi*r_p*mu*kt/dt)*g3
	);
}

//************************
void element::get_f_drag()
//************************
{
	f_drag = -6*pi*r*mu*(u-uf);
}



//************************************
void element::advect(std::string type)
//************************************
{
	
	// heun method
	if (type=="1") // first step, move particles with velocity evaluted at initial position 
	{
		u_old = uf + (f_pi + f_ext + f_pp)/6.0/pi/mu/r; 	
		p_old = p;
		p += dt*u_old; 
	}	
	else if (type=="2") //second step, evalute velocity at intermediate point, move particles using averaged velocity 
	{
		u = uf + (f_pi + f_ext + f_pp)/6.0/pi/mu/r; 	
		p = p_old + 0.5*dt*(u_old+u);	
	}
			
}


//************************ 
void element::enforce_bc()
//************************
{
	p.set_value
	(
	fmod(p.get_x()+Lx,Lx),
	fmod(p.get_y()+Ly,Ly),
	fmod(p.get_z()+Lz,Lz)
	);	
} 
