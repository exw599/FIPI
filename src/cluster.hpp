/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
 *    
 *   	https://bottogroup.wordpress.com/
 *     
 *     	Chuan Gu, c.gu@qmul.ac.uk
 *     
 *    	Copyright (2016) Botto Research Group
 *      
 * 	This software is distributed under the GNU General Public License.
 *
 * ------------------------------------------------------------------------- */


#include"head.hpp"
#include"element.hpp"
#include"linked_list.hpp"
#include"field.hpp"

#ifndef cluster_inc
#define cluster_inc

//***********
class cluster
//***********
{

        friend class list;
        friend class node;
        friend class element;

        private:
                int N;        // number of particles
                element* particle;      // array of particles
                list*    cell;          // list of particles ID for each cell
               	 
		int nx_cell;
		
		int ny_cell; 
		
		int nz_cell; 
		
		double dx_cell;
	
		double dy_cell;

		double dz_cell;

		//!
		int  it_mark;	
		double tau;

 
	public:
                cluster();

		
		~cluster();
       
		void distribute();
		
		bool check_overlap(int test_it);         
                
		void sorting();  
		
		void interpolate_local_variables
		     (const field& u, 
		      const field& c, 
                      const field& grad_c);
		
		void compute_f_pp(int t);
       	
		void compute_f_ext(int t);

		void compute_f_pi();

		void compute_f_drag();
		
		void compute_f_rand();
		         
		void advect(std::string type);

		void enforce_bc();
            
  
		const field compute_f_r();
                
		void output(std::string name);
		
		void summary();
};

#endif
