/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
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

//*********************
bool output_time(int t)
//*********************
{
        return remainder((t+1)*100,t_end)==0;
}

//**********************************************
std::string convert_to_vtk(int n, std::string name)
//**********************************************
{

        std::stringstream s;
        s<<"data/"<<name<<n<<".vtk";
        return s.str();

}

//**********************************************
std::string convert_to_csv(int n, std::string name)
//**********************************************
{
        std::stringstream s;
        s<<"data/"<<name<<n<<".csv";
        return s.str();

}

//****************************************************************************
void fftw_execute_dft(const fftw_plan& p, fftw_complex* in, fftw_complex* out);
//****************************************************************************


//*****************************
void normalize(fftw_complex* a)
//*****************************
{
        for (int it=0; it<Nx*Ny*Nz; it++)
        {
                a[it][0] /= (Nx*Ny*Nz);
                a[it][1] /= (Nx*Ny*Nz);
        }
}

//**********************************************************
void inverse(double* matrix)
//**********************************************************
{
	double a,b,c,d,e,f,g,h,i;
	double det;

	a = matrix[4]*matrix[8]-matrix[5]*matrix[7];
	b = matrix[5]*matrix[6]-matrix[3]*matrix[8]; 
	c = matrix[3]*matrix[7]-matrix[4]*matrix[6]; 
	d = matrix[2]*matrix[7]-matrix[1]*matrix[8]; 
	e = matrix[0]*matrix[8]-matrix[2]*matrix[6]; 
	f = matrix[1]*matrix[6]-matrix[0]*matrix[7]; 
	g = matrix[1]*matrix[5]-matrix[2]*matrix[4]; 
	h = matrix[2]*matrix[3]-matrix[0]*matrix[5]; 
	i = matrix[0]*matrix[4]-matrix[1]*matrix[3];

	det = matrix[0]*a+matrix[1]*b+matrix[2]*c;

	if (det == 0.0) 
	{
		std::cout<<"The matrix is singular\n";
		exit(0);
	}
	
	matrix[0] = a/det; 
	matrix[1] = d/det; 
	matrix[2] = g/det; 
	matrix[3] = b/det; 
	matrix[4] = e/det; 
	matrix[5] = h/det; 
	matrix[6] = c/det; 
	matrix[7] = f/det; 
	matrix[8] = i/det; 
}

//*****************
double sqr(double a)
//*****************
{
	return a*a;
}

//********************
double cubic(double a)
//********************
{
	return a*a*a;
}

//*******************
double sqrt(double a)
//*******************
{
	return pow(a,0.5);
}

//***************************************
double ramp(double t, double t_threshold)
//***************************************
{
	if (t < 0 ) return 0.0;
	else if( t>=0 and t<t_threshold) return t/t_threshold;
	else if( t>=t_threshold) return 1.0;
 
}


//--------------------------------------------------------------CLASS VECTOR-------------------------------------------------------------



//***********************************
vector& vector::operator+=(const vector& rhs)
//***********************************
{
	x = x + rhs.x;
	y = y + rhs.y;
	z = z + rhs.z;

	return *this;
}

//***********************************
vector& vector::operator-=(const vector& rhs)
//***********************************
{
	x = x - rhs.x;
	y = y - rhs.y;
	z = z - rhs.z;

	return *this;
}


//***********************************
vector& vector::operator=(const vector& rhs)
//***********************************
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;

	return *this;
}

//*************************************
vector& vector::operator*=(double rhs)
//*************************************
{
	x = x*rhs;
	y = y*rhs;
	z = z*rhs;
	
	return *this;
}

//*************************************
vector& vector::operator/=(double rhs)
//*************************************
{
	if (rhs == 0.0)
	{
		std::cout<<"vector/0.0\n";
		exit(0);	
	}
	
	x = x/rhs;
	y = y/rhs;
	z = z/rhs;
	
	return *this;
}

//***********************************************
const vector vector::operator+(const vector& rhs) const
//***********************************************
{
	vector output;	
	output  = *this;
	output += rhs;	
	return output;
}

//***********************************************
const vector vector::operator-(const vector& rhs) const
//***********************************************
{
	vector output;	
	output  = *this;
	output -= rhs;	
	return output;
}

//***********************************************
const double vector::operator*(const vector& rhs) const
//***********************************************
{
	double output;
	
	output = x*rhs.x + y*rhs.y + z*rhs.z;
		
	return output;
}


//**********************************************
const vector vector::operator*(double rhs) const
//**********************************************
{
	vector output;
	output  = *this;
	output *= rhs;
	return output;
}

//**********************************************
const vector vector::operator/(double rhs) const
//**********************************************
{
	vector output;
	output  = *this;
	output /= rhs;
	return output;
}

//**********************************************
const bool vector::operator==(const vector& rhs) const
//**********************************************
{
	bool output = false;
		
	if (x == rhs.x and y==rhs.y) 
	{
	if (z == rhs.z)
	{	
		output = true;
	}
	}

	return output;	
}



//-----------------------------------------------------------OUT OF CLASS FUNCTION FOR VECTOR-------------------------------------------------------


//****************************************************
const vector operator*(double lhs, const vector& rhs)
//****************************************************
{
	return rhs*lhs;
}

