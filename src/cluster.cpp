/*----------------------------------------------------------------------
 * 
 * 	FIPI - Fast Interface Particle Interactions
 *    
 *    	https://bottogroup.wordpress.com/
 *      
 *      Chuan Gu, c.gu@qmul.ac.uk
 *      
 *      Copyright (2016) Botto Research Group
 *      
 *      This software is distributed under the GNU General Public License.
 *
 * ------------------------------------------------------------------------- */

#include"head.hpp"
#include"primitive.hpp"
#include"cluster.hpp"
#include"field.hpp"

//****************
cluster::cluster()
//****************
{		
	N = N_p;		

	nx_cell = int(floor(Lx/r_c + 2*margin/r_c));
	ny_cell = int(floor(Ly/r_c + 2*margin/r_c));
	nz_cell = int(floor(Lz/r_c + 2*margin/r_c));

	dx_cell = fmod(Lx + 2*margin,r_c)/nx_cell + r_c;
	dy_cell = fmod(Ly + 2*margin,r_c)/ny_cell + r_c;
	dz_cell = fmod(Lz + 2*margin,r_c)/nz_cell + r_c;
		
	particle = new element[N];
	
	cell     = new list[nx_cell*ny_cell*nz_cell];

	tau = 0.0;
	it_mark = -1;
}


//*****************
cluster::~cluster()
//*****************
{
	delete [] particle;
	delete [] cell;
}


//*************************
void cluster::distribute()
//*************************
{ 

	double rnd1,rnd2,rnd3;

        for (int it=0; it<N; it++)
        {
                particle[it].flag = 0;

                rnd1 = double(rand()%1000)/1e3;
                rnd2 = double(rand()%1000)/1e3;
                rnd3 = double(rand()%1000)/1e3;

                particle[it].p.set_value
                (
               		Lx*rnd1,
			Ly*rnd2,
			Lz*rnd3 
		);

                while (check_overlap(it))
                {
			rnd1 = double(rand()%1000)/1e3;
			rnd2 = double(rand()%1000)/1e3;
			rnd3 = double(rand()%1000)/1e3;

			particle[it].p.set_value
			(
				Lx*rnd1,
				Ly*rnd2,
				Lz*rnd3 
			); 
		}
	
		particle[it].p_old = particle[it].p;
        }


}




//**************************************
bool cluster::check_overlap(int test_it)
//**************************************
{
	bool output = false;

	for (int it=0; it<test_it; it++)
	{
	if ((particle[it].p-particle[test_it].p).norm()<2.0*r_p)
	{
		output = true;
		break;
	}	
	}
	
	return output;
}

//**********************
void cluster::sorting()
//**********************
{
	int i,j,k;
		
	// clear all cell lists
	for (int n=0; n<nx_cell*ny_cell*nz_cell; n++)		
	{	
		cell[n].clear();	
	}


	double y_max = 0.0;

	// sorting all particles into local cells	
	for (int it=0; it<N; it++)
	{	
		//!	
		if (particle[it].p.get_y() > y_max) 
		{
			y_max = particle[it].p.get_y(); 	
			it_mark = it;								
		}
			
		i = int(floor(particle[it].p.get_x()/dx_cell + margin/dx_cell));	
		j = int(floor(particle[it].p.get_y()/dy_cell + margin/dy_cell));	
		k = int(floor(particle[it].p.get_z()/dz_cell + margin/dz_cell));	
	
		if (i < 0 or i >= nx_cell) 
		{
			std::cout<<"i="<<i<<" x out of sorting range\n";
			output("particle_error.vtk");	
			exit(0);	
		}	
		
		if (j < 0 or j >= ny_cell) 
		{
			std::cout<<"j="<<j<<" y out of sorting range\n";
			output("particle_error.vtk");	
			exit(0);
		}
		
		if (k < 0 or k >= nz_cell) 
		{			
			std::cout<<"k="<<k<<" z out of sorting range\n";
			output("particle_error.vtk");	
			exit(0);
		}
		
		cell[i+j*nx_cell+k*nx_cell*ny_cell].append(it);						
	       
	}
	
}



//*****************************************************
void cluster::compute_f_pp(int t)
//*****************************************************
{	
	int o,p,q;         // LATTICE COORDINATE OF TEST PARTICLE 
	int    co;    	   // coordinate of 27 neigboring sorting cells	
	int    id; 	   // id of particles at current node of current neighboring cell	
	node*  node_temp;  // iterator node pointer for current neigboring cell

	vector d;
	vector force_n;	

	vector bc_x(Lx,0.0,0.0);
	vector bc_y(0.0,Ly,0.0);
	vector bc_z(0.0,0.0,Lz);
		
	tau = 0.0;
	
	for (int it=0; it<N; it++) particle[it].f_pp.set_value(0.0,0.0,0.0); 		
			
	for (int it=0; it<N; it++)
	{
	

		o = int(floor(particle[it].p.get_x()/dx_cell + margin/dx_cell));
		p = int(floor(particle[it].p.get_y()/dy_cell + margin/dy_cell));
		q = int(floor(particle[it].p.get_z()/dz_cell + margin/dz_cell));
	
	
		for (int k=q-1; k<=q+1; k++)	
		{	
		for (int j=p-1; j<=p+1; j++)	
		{	
		for (int i=o-1; i<=o+1; i++)	
		{	
       	
		co  = (i+nx_cell)%nx_cell;
		co += (j+ny_cell)%ny_cell * nx_cell;
		co += (k+nz_cell)%nz_cell * nx_cell*ny_cell;
		
		node_temp = cell[co].front;	

		while(node_temp != NULL)
		{									
		
		id = node_temp->value;		

		if (id > it)
		{

			d = particle[it].p - particle[id].p;

			if      (d.get_x()> dx_cell*(nx_cell-2)) d -= bc_x;
			else if (d.get_x()<-dx_cell*(nx_cell-2)) d += bc_x;

			if      (d.get_y()> dy_cell*(ny_cell-2)) d -= bc_y;
			else if (d.get_y()<-dy_cell*(ny_cell-2)) d += bc_y;

			if      (d.get_z()> dz_cell*(nz_cell-2)) d -= bc_z;
			else if (d.get_z()<-dz_cell*(nz_cell-2)) d += bc_z;
			
			if (d.norm() == 0.0)
			{
			d.set_value
			(
			double(rand()%1000)/1e3*2.0*r_p-r_p,
			double(rand()%1000)/1e3*2.0*r_p-r_p,
			double(rand()%1000)/1e3*2.0*r_p-r_p
			);
			}

				
			//instantaneous force 	
			if (d.norm()<r_eq)  	
			{		
				force_n = (r_eq-d.norm())
					   *ramp(t,0.2*t_end)*k_n
					   *(d/d.norm());	

			particle[it].f_pp += force_n ;
			particle[id].f_pp -= force_n ;	
		

			if (particle[it].f_pi.norm() != 0.0 and particle[it_mark].f_pi.norm() != 0.0) 
			{
				tau += 0.5*force_n*d; 
			}	
			
		
			}
 
		} // if id > it

		node_temp = node_temp->next;		
		} //while node_temp != NULL  

		} //loop over all 27 adjacent cells
		} 
		}
 
	} //loop it

}


//********************************
void cluster::compute_f_ext(int t)
//********************************
{
	for (int it=0; it<N; it++) particle[it].get_f_ext(t); 
}

//**************************
void cluster::compute_f_pi()
//**************************
{
	for (int it=0; it<N; it++) particle[it].get_f_pi();		
}

//****************************
void cluster::compute_f_drag()
//****************************
{
	for (int it=0; it<N; it++) particle[it].get_f_drag();
}

//****************************
void cluster::compute_f_rand()
//****************************
{
	for (int it=0; it<N; it++) particle[it].get_f_rand();
}

//********************************************************************************************
void cluster::interpolate_local_variables(const field& u, const field& c, const field& grad_c)
//********************************************************************************************
{

	for (int it=0; it<N; it++) 
	{	
	
	if (switch_ns) particle[it].get_local_u(u);
		
	if (switch_ch)
	{ 
	particle[it].get_local_c(c);
        particle[it].get_local_dcdx(grad_c);
	}
	
	}
}

//************************************
void cluster::advect(std::string type)
//************************************
{ 
	
	for (int it=0; it<N; it++) 
	{
		particle[it].advect(type);
	}
}

//************************
void cluster::enforce_bc()
//************************
{
	for (int it=0; it<N; it++) particle[it].enforce_bc();
}

//********************************
const field cluster::compute_f_r()
//********************************
{
	vector f;

        int co;
        int x, y, z; 	//coordinates of starting node of local cell 
        double o, p, q; // distance fraction within local cell along 3 directions
        	
	field f_r("VECTORS");

   	for (int it=0; it<N; it++)
	{

	f = (particle[it].f_pp + particle[it].f_ext)/dx/dy/dz;
	
	x = int(floor(particle[it].p.get_x()/dx));
	y = int(floor(particle[it].p.get_y()/dy));
	z = int(floor(particle[it].p.get_z()/dz));

	o = particle[it].p.get_x()/dx - floor(particle[it].p.get_x()/dx);
	p = particle[it].p.get_y()/dy - floor(particle[it].p.get_y()/dy);
	q = particle[it].p.get_z()/dz - floor(particle[it].p.get_z()/dz);

	co = x + y*Nx + z*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*(1-o)*(1-p)*(1-q));
	f_r.add_value(1, co, f.get_y()*(1-o)*(1-p)*(1-q));
	f_r.add_value(2, co, f.get_z()*(1-o)*(1-p)*(1-q));

	co = (x+1)%Nx + y*Nx + z*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*   o*(1-p)*(1-q));
	f_r.add_value(1, co, f.get_y()*   o*(1-p)*(1-q));
	f_r.add_value(2, co, f.get_z()*   o*(1-p)*(1-q));

	co = (x+1)%Nx + (y+1)%Ny*Nx + z*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*   o*p*(1-q));
	f_r.add_value(1, co, f.get_y()*   o*p*(1-q));
	f_r.add_value(2, co, f.get_z()*   o*p*(1-q));

	co = (x+1)%Nx + y*Nx + (z+1)%Nz*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*   o*(1-p)*q);
	f_r.add_value(1, co, f.get_y()*   o*(1-p)*q);
	f_r.add_value(2, co, f.get_z()*   o*(1-p)*q);

	co = (x+1)%Nx + (y+1)%Ny*Nx + (z+1)%Nz*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*   o*p*q);
	f_r.add_value(1, co, f.get_y()*   o*p*q);
	f_r.add_value(2, co, f.get_z()*   o*p*q);

	co = x + y*Nx + (z+1)%Nz*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*(1-o)*(1-p)*q);
	f_r.add_value(1, co, f.get_y()*(1-o)*(1-p)*q);
	f_r.add_value(2, co, f.get_z()*(1-o)*(1-p)*q);

	co = x + (y+1)%Ny*Nx + z*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*(1-o)*p*(1-q));
	f_r.add_value(1, co, f.get_y()*(1-o)*p*(1-q));
	f_r.add_value(2, co, f.get_z()*(1-o)*p*(1-q));

	co = x + (y+1)%Ny*Nx + (z+1)%Nz*Nx*Ny;
	f_r.add_value(0, co, f.get_x()*(1-o)*p*q);
	f_r.add_value(1, co, f.get_y()*(1-o)*p*q);
	f_r.add_value(2, co, f.get_z()*(1-o)*p*q);
        }

	return f_r;
}



//************************************
void cluster::output(std::string name)
//************************************
{
        std::ofstream traject(name.c_str());
        traject<<"# vtk DataFile Version 2.0\n";
        traject<<"particles\n";
        traject<<"ASCII\n";
        traject<<"\n";

        traject<<"DATASET UNSTRUCTURED_GRID\n";
        traject<<"POINTS "<<N<<" double\n";

        for (int it=0; it<N; it++)
        {	
		traject<<particle[it].p.get_x()/dx<<" ";
		traject<<particle[it].p.get_y()/dy<<" ";
		traject<<particle[it].p.get_z()/dz<<"\n";
        }
	
	traject<<"\n";
        traject<<"POINT_DATA "<<N<<"\n";
        traject<<"SCALARS radius double\n";
        traject<<"LOOKUP_TABLE default\n";

        for (int it=0; it<N; it++)
        {
		traject<<particle[it].r/dx<<"\n";
        }

	traject<<"\n";
        traject<<"SCALARS f_pp double\n";
        traject<<"LOOKUP_TABLE default\n";

        for (int it=0; it<N; it++)
        {
		traject<<particle[it].f_pp.norm()/A/pi/sigma/r_p<<"\n";
	}
	
	traject.close();

	std::cout<<"packing energy= "<<tau<<" ";	
}


//*********************
void cluster::summary()
//************&********
{

double t_ch  = cubic(epsilon)/mobility/sigma;
double t_adv = 6*mu*r_p/A/sigma;
double u_t = 0.5*A*pi*sigma*r_p/6/pi/mu/r_p; 

std::cout<<"----------------SUMMARY------------------\n";
std::cout<<"Total simulation time = "<<dt*t_end<<"\n";
std::cout<<"Delta t = "<<dt<<"\n";
std::cout<<"\n";
std::cout<<"CH time scale = "<<t_ch<<"\n";
std::cout<<"S = "<<sqrt(mu*mobility)/epsilon<<"\n";
std::cout<<"Cn = "<<epsilon/r_drop<<"\n";
std::cout<<"l_map/r_p = "<<2*sqrt(2)*epsilon/r_p<<"\n";
std::cout<<"\n";
std::cout<<"------------------------------------------\n";
std::cout<<"Fluid viscosity mu = "<<mu<<"\n";
std::cout<<"number of particles = "<<N_p<<"\n";
std::cout<<"k_n_initial = "<<k_n<<"\n";
std::cout<<""<<"\n";
}

