/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
 *   
 *    	https://bottogroup.wordpress.com/
 *     
 *    	Chuan Gu, c.gu@qmul.ac.uk
 *     
 *    	Copyright (2016) Botto Research Group
 *     
 *	This software is distributed under the GNU General Public License.
 * 
 * ------------------------------------------------------------------------- */



#include"head.hpp"
#include"primitive.hpp"
#include"field.hpp"

//----------------------------------------------------CONSTRUCTORS AND DESTRUCTORS-----------------------------------------

//********************************
field::field(std::string type_id)
//********************************
{
	type = type_id;					
	nx = Nx;
	ny = Ny;
	nz = Nz;
	
	if      (type == "SCALARS") dimension = 1; 
	else if (type == "VECTORS") dimension = 3; 
	else if (type == "TENSORS") dimension = 9;
	else
	{
		std::cout<<"UNDEFINED FIELD TYPE\n";
		exit(0);
	}
	
	v =  new fftw_complex*[dimension];

	for (int d=0; d<dimension; d++) 
	{
		v[d] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
	}
	
	for (int d=0; d<dimension; d++)
	{
	for (int it=0; it<nx*ny*nz; it++) 
	{	
		v[d][it][0] = 0.0;
		v[d][it][1] = 0.0;
	}	
	}
	
	in_fspace = false;
}


//****************************************************
field::field(std::string type_id, std::string name_id)
//****************************************************
{
	type = type_id;					
	nx = Nx;
	ny = Ny;
	nz = Nz;

	if      (type == "SCALARS") dimension = 1; 
	else if (type == "VECTORS") dimension = 3; 
	else if (type == "TENSORS") dimension = 9;

	v =  new fftw_complex*[dimension];

	for (int d=0; d<dimension; d++)
	{
		v[d] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
	}
	
	in_fspace = false;

	int i,j,k;

	if (name_id=="phase_field")
	{	
		//************** initialise phase field here ******************************	
		for (int it=0; it<nx*ny*nz; it++) 
		{
			k = it/(nx*ny);
			j = it%(nx*ny)/nx;
			i = it - nx*j - nx*ny*k;		  		
		//	v[0][it][0] = tanh((sqrt(sqr(k*dz-Lz/2.0))-0.3*Lx)/sqrt(2)/epsilon);
		//	v[0][it][0] = tanh((k*dz-Lz/2.0)/sqrt(2)/epsilon);
	        //      v[0][it][0] = tanh((sqrt(sqr(i*Lx/Nx-0.5*Lx)+sqr(j*Ly/Ny-0.5*Ly)+sqr(k*Lz/Nz-0.5*Lz))-r_drop)/sqrt(2.0)/epsilon);	
			v[0][it][0] = double(rand()%1000)/1e3*0.2 - 0.1;	
			v[0][it][1] = 0.0;
		}	
	}	
	else if (name_id=="wave_vector")
	{	
		for (int it=0; it<nx*ny*nz; it++) 
		{
			k = it/(nx*ny);
			j = it%(nx*ny)/nx;
			i = it%(nx*ny)%nx;		  	
			
			if      (i <=  nx/2) v[0][it][0] = i*2.0*pi/Lx;
			else                 v[0][it][0] = (i - nx)*2.0*pi/Lx;

			if      (j <=  ny/2) v[1][it][0] = j*2.0*pi/Ly;
			else                 v[1][it][0] = (j - ny)*2.0*pi/Ly;

			if      (k <=  nz/2) v[2][it][0] = k*2.0*pi/Lz;
			else                 v[2][it][0] = (k - nz)*2.0*pi/Lz;
			
			v[0][it][1]  = 0.0;
			v[1][it][1]  = 0.0;
			v[2][it][1]  = 0.0;		
		}
	}
	else if(name_id=="gravity")
	{	
		for (int it=0; it<nx*ny*nz; it++)
		{
			v[0][it][0] = 0.0;	
			v[0][it][1] = 0.0;	
			v[1][it][0] = 0.0;	
			v[1][it][1] = 0.0;	
			v[2][it][0] = gravity;	
			v[2][it][1] = 0.0;	
		} 
	}	
	else
	{
		std::cout<<"UNDEFINED FIELD NAME\n";
		exit(0);	
	}
}

//*****************************************
void field::set_value(int d,int it, double value)
//*****************************************
{
	v[d][it][0] = value;
	v[d][it][1] = 0.0;
}

//*****************************************
void field::add_value(int d,int it, double value)
//*****************************************
{
	v[d][it][0] += value;
	v[d][it][1] = 0.0;
}


//*************
field::~field()
//*************
{
	free();	
	delete [] v;
}

//**********************************
void field::malloc(std::string space)
//**********************************
{
        for (int d=0; d<dimension; d++) v[d] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);

	if      (space == "f") in_fspace = true;
	else if (space == "p") in_fspace = false;
}



//**********************************
void field::reset(std::string space)
//**********************************
{	
	for (int d=0; d<dimension; d++)
	{
	for (int it=0; it<nx*ny*nz; it++) 
	{	
		v[d][it][0] = 0.0;
		v[d][it][1] = 0.0;
	}	
	}

	if      (space == "f") in_fspace = true;
	else if (space == "p") in_fspace = false;
}

//*******************
const field field::norm() const
//*******************
{

        if (in_fspace)
        {
                std::cout<<"COMPUTING NORM IN FREQUENCY SPACE\n";
                exit(0);
        }

	field output("SCALARS");
		
	for (int it=0; it<nx*ny*nz; it++) 
	{	
		for(int d=0; d<dimension; d++)
		{	
			output.v[0][it][0] += sqr(v[d][it][0]);	
		}
		
		output.v[0][it][0] = sqrt(output.v[0][it][0]);	
		output.v[0][it][1] = 0.0;						
	}

	return output;
}

//*************************************************
double field::interpolate(const vector& position, int d) const
//*************************************************
{

        double value_off;

        int i = int(floor(position.get_x()/dx));
        int j = int(floor(position.get_y()/dy));
        int k = int(floor(position.get_z()/dz));

        double m = position.get_x()/dx - floor(position.get_x()/dx);
        double n = position.get_y()/dy - floor(position.get_y()/dy);
        double o = position.get_z()/dz - floor(position.get_z()/dz);

        value_off  = (1.0-m)*(1.0-n)*(1.0-o)*v[d][ i       +  j      *Nx + k*Nx*Ny][0];
        value_off += (1.0-m)*     n *(1.0-o)*v[d][ i       + (j+1)%Ny*Nx + k*Nx*Ny][0];
        value_off +=      m *(1.0-n)*(1.0-o)*v[d][(i+1)%Nx +  j      *Nx + k*Nx*Ny][0];
        value_off +=      m *     n *(1.0-o)*v[d][(i+1)%Nx + (j+1)%Ny*Nx + k*Nx*Ny][0];

        value_off += (1.0-m)*(1.0-n)*     o *v[d][ i       +  j      *Nx + (k+1)%Nz*Nx*Ny][0];
        value_off += (1.0-m)*     n *     o *v[d][ i       + (j+1)%Ny*Nx + (k+1)%Nz*Nx*Ny][0];
        value_off +=      m *(1.0-n)*     o *v[d][(i+1)%Nx +  j      *Nx + (k+1)%Nz*Nx*Ny][0];
        value_off +=      m *     n *     o *v[d][(i+1)%Nx + (j+1)%Ny*Nx + (k+1)%Nz*Nx*Ny][0];

        return value_off;
}

//*************************************************
double field::interpolate2(const vector& position, int d) const
//*************************************************
{

        double value_off;
        int size = 3; //interpolation kernel size in grid cells 

        int i = int(floor(position.get_x()/dx));
        int j = int(floor(position.get_y()/dy));
        int k = int(floor(position.get_z()/dz));

        double m = position.get_x()/dx - floor(position.get_x()/dx);
        double n = position.get_y()/dy - floor(position.get_y()/dy);
        double o = position.get_z()/dz - floor(position.get_z()/dz);

        m = m/size + 1.0/size;
        n = n/size + 1.0/size;
        o = o/size + 1.0/size;

        value_off  = (1.0-m)*(1.0-n)*(1.0-o)*v[d][ (i-(1-size)/2 + Nx)%Nx + (j-(1-size)/2 + Ny)%Ny*Nx + (k-(1-size)/2 + Nz)%Nz*Nx*Ny][0];
        value_off += (1.0-m)*     n *(1.0-o)*v[d][ (i-(1-size)/2 + Nx)%Nx + (j+(1+size)/2     )%Ny*Nx + (k-(1-size)/2 + Nz)%Nz*Nx*Ny][0];
        value_off +=      m *(1.0-n)*(1.0-o)*v[d][ (i+(1+size)/2     )%Nx + (j-(1-size)/2 + Ny)%Ny*Nx + (k-(1-size)/2 + Nz)%Nz*Nx*Ny][0];
        value_off +=      m *     n *(1.0-o)*v[d][ (i+(1+size)/2     )%Nx + (j+(1+size)/2     )%Ny*Nx + (k-(1-size)/2 + Nz)%Nz*Nx*Ny][0];

        value_off += (1.0-m)*(1.0-n)*     o *v[d][ (i-(1-size)/2 + Nx)%Nx + (j-(1-size)/2 + Ny)%Ny*Nx + (k+(1+size)/2)%Nz*Nx*Ny][0];
        value_off += (1.0-m)*     n *     o *v[d][ (i-(1-size)/2 + Nx)%Nx + (j+(1+size)/2     )%Ny*Nx + (k+(1+size)/2)%Nz*Nx*Ny][0];
        value_off +=      m *(1.0-n)*     o *v[d][ (i+(1+size)/2     )%Nx + (j-(1-size)/2 + Ny)%Ny*Nx + (k+(1+size)/2)%Nz*Nx*Ny][0];
        value_off +=      m *     n *     o *v[d][ (i+(1+size)/2     )%Nx + (j+(1+size)/2     )%Ny*Nx + (k+(1+size)/2)%Nz*Nx*Ny][0];

        return value_off;
}



//-----------------------------------------------------------------VECTOR CALCULUS------------------------------------------------------

//********************************************
const field field::grad(const field& kf) const
//********************************************
{
	if (not in_fspace)
        {
                std::cout<<"YOU CAN ONLY CALCULATE GRADIENT IN FREQUENCY SPACE!\n";
                exit(0);
        }

	std::string type_of_gradient;
	
	if      (type=="SCALARS") type_of_gradient = "VECTORS";
	else if (type=="VECTORS") type_of_gradient = "TENSORS";
	
	field output(type_of_gradient); 
	
	output.in_fspace = true;
			
	for (int d=0; d<output.dimension; d++)		
	{
	for (int it=0; it<nx*ny*nz; it++)
	{ 
		if ( 	
			(it%(nx*ny)%nx == nx/2 and d%3 == 0) 
		     or	(it%(nx*ny)/nx == ny/2 and d%3 == 1)
		     or	(it/(nx*ny)    == nz/2 and d%3 == 2)
		)
		{
			output.v[d][it][0] = 0.0;
			output.v[d][it][1] = 0.0;
		}
		else
		{
			output.v[d][it][0] = -kf.v[d%3][it][0]*v[d/3][it][1];
			output.v[d][it][1] =  kf.v[d%3][it][0]*v[d/3][it][0];
		}
	}
	}	
	
	return output;
}


//*******************************************
const field field::div(const field& kf) const
//*******************************************
{
        if (not in_fspace)
        {
                std::cout<<"YOU CAN ONLY CALCULATE DIVERGENCE IN FREQUENCY SPACE!\n";
                exit(0);
        }
	
	std::string type_of_divergence;

	if      (type=="SCALARS") type_of_divergence = "SCALARS"; 
	else if (type=="VECTORS") type_of_divergence = "SCALARS";
        else if (type=="TENSORS") type_of_divergence = "VECTORS";
 	
	field output(type_of_divergence); 
			
	output.in_fspace = true;

        for (int d=0; d<dimension; d++)
        {
        for (int it=0; it<nx*ny*nz; it++)
	{ 
		output.v[d/3][it][0] += - kf.v[d%3][it][0] * v[d][it][1];
		output.v[d/3][it][1] +=   kf.v[d%3][it][0] * v[d][it][0];
	} 
	}
	
	return output;
}

//*******************************************
const field field::dot(const field& b) const
//*******************************************
{	
 	if (dimension != b.dimension)
        {
                std::cout<<"dot() TWO FIELDS WITH DIFFERENT DIMENSIONS";
                exit(0);
        }

        field output("SCALARS"); 
	
	if (in_fspace or b.in_fspace ) 
	{		
		output.in_fspace = true;

		for (int it=0; it<nx*ny*nz; it++)
		{
		for (int d=0; d<dimension; d++)
		{	
			output.v[0][it][0]  += v[d][it][0]*b.v[d][it][0] - v[d][it][1]*b.v[d][it][1];
			output.v[0][it][1]  += v[d][it][1]*b.v[d][it][0] + v[d][it][0]*b.v[d][it][1];
		}
        	}
	}	
	else 
	{
		output.in_fspace = false;
		
		for (int it=0; it<nx*ny*nz; it++)
		{
		for (int d=0; d<dimension; d++)
		{	
			output.v[0][it][0]  += v[d][it][0]*b.v[d][it][0];

		}
		output.v[0][it][1] = 0.0;
        	}
	}
        
	return output;

}

//***************************************
const field field::power(double b) const
//***************************************
{
	if (in_fspace) 
	{
		std::cout<<"POWER() FIELD IN FREQUENCY SPACE\n";
		exit(0);
	}	

        field output("SCALARS"); 
		
	output.in_fspace = false;

       	for (int d=0; d<dimension; d++)
      	{
        for (int it=0; it<output.nx*output.ny*output.nz; it++)
        {
		output.v[0][it][0] += pow(v[d][it][0],b);
		output.v[0][it][1]  = 0.0; 
	}		
        }
        
	return output;
}



//------------------------------------------------------------------FFT------------------------------------------------------------

//*****************************************
void field::fft(const fftw_plan& p_forward)
//*****************************************
{
	for (int d=0; d<dimension; d++) fftw_execute_dft(p_forward, v[d], v[d]);
	in_fspace = true;
}

//*******************************************
void field::ifft(const fftw_plan& p_backward)
//*******************************************
{
	for (int d=0; d<dimension; d++) 
	{		
		fftw_execute_dft(p_backward, v[d], v[d]);
		normalize(v[d]);
	}	
	in_fspace = false;
}

//****************
void field::free()
//****************
{
 	for (int d=0; d<dimension; d++) fftw_free(v[d]);
}



//--------------------------------------------------------------OUTPUT--------------------------------------------------------------------------

//**************************************************************************
void field::output_to_file(std::string file_name, std::string variable_name)
//**************************************************************************
{

	if (in_fspace)
	{
		std::cout<<"YOU ARE WRITING FIELD IN FREQUENCY SPACE INTO FILES\n";
		exit(0);
	}

        std::ofstream output(file_name.c_str());
        output<<"# vtk DataFile Version 2.0\n";
        output<<variable_name<<" PROFILE\n";
        output<<"ASCII\n";
        output<<"\n";
        output<<"DATASET STRUCTURED_POINTS\n";
        output<<"DIMENSIONS "<<nx<<" "<<ny<<" "<<nz<<"\n";
        output<<"ASPECT_RATIO 1 1 1\n";
        output<<"ORIGIN 0 0 0\n";
        output<<"POINT_DATA "<<nx*ny*nz<<"\n";
        output<<type<<" "<<variable_name<<" double\n";
        
	if (dimension==1) output<<"LOOKUP_TABLE default\n";


	for (int it=0; it<nx*ny*nz; it++)
	{ 
	for (int d=0; d<dimension; d++) 
        { 
		if (d==dimension-1) output<<v[d][it][0]<<"\n";
		else                output<<v[d][it][0]<<" ";
	}	
        }

        output.close();


	//!
	double sum = 0.0; 
	if (variable_name == "free_energy")
	{
		for (int it=0; it<nx*ny*nz; it++)
		{
			sum += v[0][it][0]*cubic(dx);	
		}
		
		std::cout<<"total free energy = "<<sum<<" ";	
	}

}





//----------------------------------------------------COMPOUND ASSIGNMENT-----------------------------------------



//****************************************
field& field::operator+=(const field& rhs)
//****************************************
{
   	if (dimension != rhs.dimension)
        {
                std::cout<<"+=FIELD TWO FIELDS WITH DIFFERENT DIMENSIONS\n";
                exit(0);
        }
   	
	if (in_fspace != rhs.in_fspace)
        {
                std::cout<<"+=FIELD TWO FIELDS IN DIFFERENT SPACE\n";
                exit(0);
        }


	if (in_fspace)
	{
		for (int it=0; it<nx*ny*nz; it++)
		{
		for (int d=0; d<dimension; d++)
		{

			v[d][it][0] += rhs.v[d][it][0];
			v[d][it][1] += rhs.v[d][it][1];
		}
        	}
	}
	else
	{
		for (int it=0; it<nx*ny*nz; it++)
		{
		for (int d=0; d<dimension; d++)
		{
			v[d][it][0] += rhs.v[d][it][0];
			v[d][it][1] = 0.0;	
		}	
		}

	}


        return *this;
}


//**********************************
field& field::operator+=(double rhs)
//**********************************
{
	if (in_fspace)
	{
		std::cout<<"+=DOUBLE ADDING REAL NUMBER IN FREQUENCY SPACE\n";
		exit(0);
	}
	
	for (int d=0; d<dimension; d++)
        {
        for (int it=0; it<nx*ny*nz; it++)
	{
		v[d][it][0] += rhs;
		v[d][it][1] = 0.0;
	}
        }

	return *this;
}

//****************************************
field& field::operator-=(const field& rhs)
//****************************************
{
  	if (dimension != rhs.dimension)
        {
                std::cout<<"-=FIELD TWO FIELDS WITH DIFFERENT DIMENSIONS\n";
                exit(0);
        }

	if (in_fspace != rhs.in_fspace)
        {
                std::cout<<"-=FIELD TWO FIELDS IN DIFFERENT SPACE\n";
                exit(0);
        }

	if (in_fspace)
	{	
		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			v[d][it][0] -= rhs.v[d][it][0];
			v[d][it][1] -= rhs.v[d][it][1];
		}
		}
	}
	else
	{
		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			v[d][it][0] -= rhs.v[d][it][0];
			v[d][it][1]  = 0.0;
		}
		}

	}

        return *this;
}

//**********************************
field& field::operator-=(double rhs)
//**********************************
{
	if (in_fspace)
	{
		std::cout<<"-=DOUBLE MINUS REAL NUMBER IN FREQUENCY SPACE\n";
		exit(0);
	}

	for (int d=0; d<dimension; d++)
        {
        for (int it=0; it<nx*ny*nz; it++)
	{
		v[d][it][0] -= rhs;
		v[d][it][1]  = 0.0;
	}
        } 
	return *this;
}




//****************************************
field& field::operator*=(const field& rhs)
//****************************************
{

	double cross_real,cross_imagine;

	if(rhs.dimension == 1)
	{	
		if(in_fspace or rhs.in_fspace)		
		{	
                	for (int d=0; d<dimension; d++)
        		{
        		for (int it=0; it<nx*ny*nz; it++)
                	{
				cross_real    = -v[d][it][1]*rhs.v[0][it][1];	 
				cross_imagine =  v[d][it][0]*rhs.v[0][it][1];	
				
				v[d][it][0] *= rhs.v[0][it][0];
				v[d][it][0] += cross_real;

				v[d][it][1] *= rhs.v[0][it][0];
				v[d][it][1] += cross_imagine; 
                	}
        		}
		}
		else
                {	
			for (int d=0; d<dimension; d++)
        		{
        		for (int it=0; it<nx*ny*nz; it++)
                	{	
				v[d][it][0] *= rhs.v[0][it][0];
				v[d][it][1]  = 0.0;	
			}
			}
		}	
	}
	else
	{
		std::cout<<"*= RHS FIELD IS NOT A SCALAR ONE\n";
		exit(0);	
	}
        
	return *this;
}

//**********************************
field& field::operator*=(double rhs)
//**********************************
{
	if (in_fspace)
	{	
		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			v[d][it][0]  *= rhs;
			v[d][it][1]  *= rhs;
		}
		}
	}
	else
	{
		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			v[d][it][0]  *= rhs;
			v[d][it][1]   = 0.0;
		}
		}	
	}
		
	return *this;
}

//****************************************
field& field::operator/=(const field& rhs)
//****************************************
{
  	if (rhs.in_fspace)
        {
                std::cout<<"/= RHS FIELD IN FREQUENCY SPACE\n";
                exit(0);
        }

	if (rhs.dimension ==1)
	{         
		if (in_fspace)
		{	 
			for (int d=0; d<dimension; d++)
        		{
        		for (int it=0; it<nx*ny*nz; it++)
                	{
                        	if (rhs.v[0][it][0] != 0.0)
                        	{
                                	v[d][it][0] /= rhs.v[0][it][0];
                                	v[d][it][1] /= rhs.v[0][it][0];
                        	}
               			else if (rhs.v[0][it][0] == 0.0)
				{
					v[d][it][0] = 0.0;	
					v[d][it][1] = 0.0;	
				} 	
			}
        		}
		}
		else
		{
			for (int d=0; d<dimension; d++)
        		{
        		for (int it=0; it<nx*ny*nz; it++)
                	{
                        	if      (rhs.v[0][it][0] != 0.0) v[d][it][0] /= rhs.v[0][it][0]; 
               			else if (rhs.v[0][it][0] == 0.0) v[d][it][0]  = 0.0;	
	
				v[d][it][1] = 0.0;	
			}
        		}
		}
	}
	else
	{
		std::cout<<"/= RHS FIELD IS NOT A SCALAR ONE\n";
		exit(0);	
	}

        return *this;
}

//**********************************
field& field::operator/=(double rhs)
//**********************************
{
	if (in_fspace)
	{
		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			if(rhs!= 0.0)
			{               
				v[d][it][0] /= rhs;
				v[d][it][1] /= rhs;
			}
			else if(rhs == 0.0)
			{
				v[d][it][0] = 0.0;
				v[d][it][1] = 0.0;
			} 
		}
		}
	}
	else
	{
		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			if     (rhs != 0.0) v[d][it][0] /= rhs;
			else if(rhs == 0.0) v[d][it][0]  = 0.0;
				
			v[d][it][1] = 0.0;
		}
		}
	}
	
	return *this;
}

//***************************************
field& field::operator=(const field& rhs)
//***************************************
{

	if (this == &rhs) 
	{		
		return *this;
	}	
	else
	{	
		if (dimension != rhs.dimension)
		{
			std::cout<<"= RHS AND LHS FIELD NOT OF THE SAME TYPE\n";	
			exit(0);	
		}

		in_fspace = rhs.in_fspace;

		for (int d=0; d<dimension; d++)
		{
		for (int it=0; it<nx*ny*nz; it++)
		{
			v[d][it][0] = rhs.v[d][it][0]; 
			v[d][it][1] = rhs.v[d][it][1]; 
		}
		}
		
		return *this;
	}	

}

//***********************
void field::print(int it)
//***********************
{
	for (int d=0; d<dimension; d++) std::cout<<v[d][it][0]<<" ";
	std::cout<<"\n";
}
//----------------------------------------------------BINARY ARITHMETIC OPERATORS-----------------------------------------





//**************************************************
const field field::operator+(const field& rhs) const
//**************************************************
{
	field result(type);
	result  = *this;
	result += rhs;
	return result;

}

//*******************************************
const field field::operator+(double rhs) const
//*******************************************
{
	field result(type);
	result  = *this; 
	result += rhs;
	return result;
}


//************************************************
const field field::operator-(const field& rhs) const
//************************************************
{
	field result(type);
	result  = *this;
	result -= rhs;
	return result;
}

//*******************************************
const field field::operator-(double rhs) const
//*******************************************
{
	field result(type);
	result  = *this;
	result -= rhs;
	return result;
}

//**************************************************
const field field::operator*(const field& rhs) const
//**************************************************
{  
        
	std::string type_of_output;
        
	if      (dimension*rhs.dimension == 1) type_of_output = "SCALARS";
        else if (dimension*rhs.dimension == 3) type_of_output = "VECTORS";
        else if (dimension*rhs.dimension == 9) type_of_output = "TENSORS";

        field output(type_of_output); 

	if (in_fspace or rhs.in_fspace) output.in_fspace = true;
	else                            output.in_fspace = false;

	if (output.in_fspace)
	{
		for (int it=0; it<nx*ny*nz; it++)
		{
		for (int d1=0; d1<dimension; d1++)
		{
		for (int d2=0; d2<rhs.dimension; d2++)
		{
			output.v[(d1+1)*(d2+1)-1][it][0] = v[d1][it][0]*rhs.v[d2][it][0] - v[d1][it][1]*rhs.v[d2][it][1];
			output.v[(d1+1)*(d2+1)-1][it][1] = v[d1][it][1]*rhs.v[d2][it][0] + v[d1][it][0]*rhs.v[d2][it][1];	
		}
		}
		}
	}
	else
	{
		for (int it=0; it<nx*ny*nz; it++)
		{
		for (int d1=0; d1<dimension; d1++)
		{
		for (int d2=0; d2<rhs.dimension; d2++)
		{
			output.v[(d1+1)*(d2+1)-1][it][0] = v[d1][it][0]*rhs.v[d2][it][0];
			output.v[(d1+1)*(d2+1)-1][it][1] = 0.0;		
		}
		}
		}
	}

        return output;
}


//******************************************
const field field::operator*(double rhs) const 
//******************************************
{
	field result(type);
	result  = *this;
	result *= rhs;
	return result;
}


//*************************************************
const field field::operator/(const field& rhs) const
//*************************************************
{
	if (rhs.in_fspace)
	{
		std::cout<<"/FIELD DIVISOR FIELD IN FREQUENCY SPACE\n";
		exit(0);
	}

	std::string type_of_output;
	if      (dimension*rhs.dimension == 1) type_of_output = "SCALARS";
        else if (dimension*rhs.dimension == 3) type_of_output = "VECTORS";
        else if (dimension*rhs.dimension == 9) type_of_output = "TENSORS";

        field output(type_of_output); 
	
	output.in_fspace = in_fspace;
               
	if (output.in_fspace)
	{ 
		for (int d1=0; d1<dimension; d1++)
		{		
		for (int d2=0; d2<rhs.dimension; d2++)	
		{	
		for (int it=0; it<nx*ny*nz; it++)
		{
			if (rhs.v[d2][it][0] != 0.0)
			{			
				output.v[(d1+1)*(d2+1)-1][it][0] = v[d1][it][0]/rhs.v[d2][it][0];
				output.v[(d1+1)*(d2+1)-1][it][1] = v[d1][it][1]/rhs.v[d2][it][0];
			}	
			else if (rhs.v[d2][it][0] == 0.0)
			{
				output.v[(d1+1)*(d2+1)-1][it][0] = 0.0;
				output.v[(d1+1)*(d2+1)-1][it][1] = 0.0;

			}
		}	
		} 
        	}
	}
	else
	{
		for (int d1=0; d1<dimension; d1++)
		{		
		for (int d2=0; d2<rhs.dimension; d2++)	
		{	
		for (int it=0; it<nx*ny*nz; it++)
		{
		
			if      (rhs.v[d2][it][0] != 0.0) output.v[(d1+1)*(d2+1)-1][it][0] = v[d1][it][0]/rhs.v[d2][it][0];
			else if (rhs.v[d2][it][0] == 0.0) output.v[(d1+1)*(d2+1)-1][it][0] = 0.0;
				
			output.v[(d1+1)*(d2+1)-1][it][1] = 0.0;
		}
		}
		}
	}
        
	return output;
}

//*******************************************
const field field::operator/(double rhs) const 
//*******************************************
{
	field result(type);
	result  = *this;	
	result /= rhs;
	return result;
}


//------------------------------------------------NON MEMBER FUNCTIONS---------------------------------



//*********************************************
const field operator+(double lhs, const field& rhs)
//*********************************************
{
	return rhs+lhs;
}

//*********************************************
const field operator*(double lhs, const field& rhs)
//*********************************************
{
	return rhs*lhs;
}

//*********************************************
const field operator-(double lhs, const field& rhs)
//*********************************************
{
	return rhs*(-1.0) + lhs;
}
