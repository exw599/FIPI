/* ----------------------------------------------------------------------
 * 	FIPI - Fast Interface Particle Interactions
 *    
 *    	https://bottogroup.wordpress.com/
 *      
 *      Chuan Gu, c.gu@qmul.ac.uk
 *      
 *      Copyright (2016) Botto Research Group
 *      
 *      This software is distributed under the GNU General Public License.	
 *------------------------------------------------------------------------- */

#include"head.hpp"
#include"field.hpp"
#include"primitive.hpp"
#include"cluster.hpp"

//***********
int main()
//***********
{

//------------------------------
srand(time(NULL)); 
clock_t start, finish; 
start = clock();

//-----------wave vector---------
field kf("VECTORS","wave_vector"); 
field kf2("SCALARS");
kf2 = kf.power(2.0);

//-------------------------------
field c("SCALARS","phase_field");
field c_new("SCALARS");
field c_old("SCALARS");
field f("SCALARS");
field f_old("SCALARS");
field grad_c("VECTORS");
field conv_c("SCALARS");
field conv_c_old("SCALARS");

c.output_to_file(convert_to_vtk(0,"phase"),"phase");
//!
field free_energy("SCALARS");


//-------------------------------
field u("VECTORS");
field pch("SCALARS"); // for surface tension 
field g("VECTORS","gravity"); // for gravity
field f_total("VECTORS");

u.output_to_file(convert_to_vtk(0,"u"),"u");


//-------------------------------
cluster suspension;
suspension.summary();

if (switch_particle) 
{
	suspension.distribute();
	suspension.output(convert_to_vtk(0,"particle"));
	std::cout<<"\n";
}


//-------------------------------
fftw_init_threads(); fftw_plan_with_nthreads(4);
fftw_complex* dummy  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
fftw_plan p_forward  = fftw_plan_dft_3d(Nz, Ny, Nx, dummy, dummy, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan p_backward = fftw_plan_dft_3d(Nz, Ny, Nx, dummy, dummy, FFTW_BACKWARD,FFTW_ESTIMATE);

for (int t=0; t<t_end; t++ ) 
{ 

	// **************************** Phase Field ***************************	
	if (switch_ch)
	{			
		f  = (lambda/sqr(epsilon))*(c.power(3.0)-c);

		f.fft(p_forward);
		c.fft(p_forward);
			
		if (switch_ns)
		{
			grad_c = c.grad(kf);	
			grad_c.ifft(p_backward);			
				
			conv_c = u.dot(grad_c);								
			conv_c.fft(p_forward);	
		}	

		if (t == 0) 
		{

			c_old = c;
			f_old = f;
			
			if (switch_ns) 
			{
				conv_c_old = conv_c;
			}
	
			c -= dt*mobility*kf2*f;					
			
			if (switch_ns) 
			{
				c -= dt*conv_c;
			}
	
			c /= 1.0 + (dt*mobility*lambda) * (kf2.power(2.0));
		} 
		else
		{ 
			c_new  = 4.0*c - c_old;	
			c_new -= (4.0*dt*mobility) * kf2 * f;
			c_new += (2.0*dt*mobility) * kf2 * f_old;
		
			if (switch_ns) 
			{
				c_new -= (4.0*dt)*conv_c;
				c_new += (2.0*dt)*conv_c_old;		
			}
			
			c_new /= 3.0 + (2.0*dt*mobility*lambda) * (kf2.power(2.0));		
			
			c_old = c;
			f_old = f; 
	
			if (switch_ns) 
			{
				conv_c_old = conv_c;	
			}
	
			c = c_new;
		}
		
		c.ifft(p_backward);	
	}



	// ****************************** Fluid Velocity Field ************************************       	
	if (switch_ns)
	{

		f_total.reset("p");
		
		if (switch_ch)
		{			
			f  = (lambda/sqr(epsilon))*(c.power(3.0)-c);	
			f.fft(p_forward);
			c.fft(p_forward);

			pch = lambda*kf2*c + f;	
			pch.ifft(p_backward);

			grad_c = c.grad(kf);
			grad_c.ifft(p_backward);	
			c.ifft(p_backward);
 	
			f_total += pch*grad_c;
		}

		if (switch_particle)
		{				
			f_total += suspension.compute_f_r();	
		}			 	
	
		f_total += ((1.0+c)/2.0*rho1 + (1.0-c)/2.0*rho2)*g;		
	
		//---------------------------------	
		f_total.fft(p_forward);		
	
		u = f_total - kf*(kf.dot(f_total))/kf2;
		u /= mu*kf2;
		u.ifft(p_backward);		

	}


	
	//********************** Particles **********************************	
	if (switch_particle)
	{	
		if (switch_ch)
		{
			c.fft(p_forward);
			grad_c = c.grad(kf);	
			grad_c.ifft(p_backward);	
			c.ifft(p_backward);	
		}	
		
		//------------------------------------	

		suspension.interpolate_local_variables(u, c, grad_c);		
		suspension.sorting();

		suspension.compute_f_pi();
		suspension.compute_f_pp(t+1);
		suspension.compute_f_ext(t+1);
	
		suspension.advect("1");			
		suspension.enforce_bc();							
	
	
		suspension.interpolate_local_variables(u, c, grad_c);						
		suspension.sorting();

		suspension.compute_f_pi();
		suspension.compute_f_pp(t+1);
		suspension.compute_f_ext(t+1);
		
		suspension.advect("2");		
		suspension.enforce_bc();							

	}


	
	
	// ****************** OUTPUT AT EVERY PERCENTAGE ******************* 
	if (output_time(t)) 
	{
		std::cout<<"Time = "<<(t+1)*dt<<" ";
	
		if (switch_ch)
		{
			c.output_to_file(convert_to_vtk(t+1,"phase"),"phase");
			free_energy = lambda*(0.5*grad_c.power(2.0) + 0.25/sqr(epsilon)*(c.power(2.0)-1.0).power(2.0));	
			free_energy.output_to_file(convert_to_vtk(t+1,"free_energy"),"free_energy");	
		}
	
	//	if (switch_ns) u.output_to_file(convert_to_vtk(t+1,"u"),"u");	
		
		if (switch_particle) suspension.output(convert_to_vtk(t+1,"particle"));			
	
		std::cout<<"\n";
	} 

}//time loop

fftw_destroy_plan(p_forward);
fftw_destroy_plan(p_backward);

finish = clock();
std::cout<<"Wall clock time: "<<(finish-start)/(CLOCKS_PER_SEC*3600)<<" Hour "<<(finish-start)%(CLOCKS_PER_SEC*3600)/(60*CLOCKS_PER_SEC)<<" Min\n";

return 0;
}
