/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
 *   
 *    	https://bottogroup.wordpress.com/
 *     
 *      Chuan Gu, c.gu@qmul.ac.uk
 *     
 *     	Copyright (2016) Botto Research Group
 *  
 *	This software is distributed under the GNU General Public License.
 *  
 * ------------------------------------------------------------------------- */

#include"field.hpp"
#include"primitive.hpp"
#ifndef element_inc
#define element_inc

//***********
class element
//***********
{
	friend class cluster;

        private:
 
		double r;   //radius
                
		double rho; //density                     
                
		vector p;   //position
 		vector p_old;
               
		vector u;   //velocity
		vector u_old;                
                
		int  flag;     

		//************************************

                vector f_pp;  //total particle-particle
                
		vector f_rand;//random force
                
		vector f_pi;  //particle-interface
                
		vector f_drag;//hydrodynamic drag
                
		vector f_ext; //external force
		
		//************************************
 
                double c;     //local phase field
                
		vector dcdx;  //local phase field gradient
                
		vector uf;    //local fluid velocity

        public:
                element();

                double volume(){return 4.0/3.0*pi*pow(r,3);}
		
		//***************
                
                void get_local_u(const field& u);
                
		void get_local_c(const field& c);
                
		void get_local_dcdx(const field& grad_c);


		//***************

                void get_f_rand();
                
		void get_f_pi();
                
		void get_f_drag();
                
		void get_f_ext(int t);
		

		
		//***************
	
                void advect(std::string type);
		
		void enforce_bc();
		
		void wall_interaction(vector v_wall, vector norm_wall);
};
#endif
