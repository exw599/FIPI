/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
 *    
 *    	https://bottogroup.wordpress.com/
 *     
 *    	Chuan Gu, c.gu@qmul.ac.uk
 *     
 *   	Copyright (2016) Botto Research Group
 *  
 *	This software is distributed under the GNU General Public License.
 * 
 * ------------------------------------------------------------------------- */

#ifndef primitive_inc 
#define primitive_inc 

bool output_time(int t);

std::string convert_to_vtk(int n, std::string name);

std::string convert_to_csv(int n, std::string name);

void normalize(fftw_complex* a);

void inverse(double* matrix);

double sqr(double a);

double cubic(double a);

double sqrt(double a);

double ramp(double t, double t_threshold);

//**********
class vector
//**********
{
private:
	double x;
	double y;
	double z;
public:
	vector()
	: x(0.0), y(0.0), z(0.0)
	{}

	vector(double a, double b, double c)
	: x(a), y(b), z(c)
	{}

	void set_value(double a, double b, double c){x=a; y=b; z=c;} 
	
	double get_x() const {return x;}
	double get_y() const {return y;}
	double get_z() const {return z;}
	double norm() const  {return pow(pow(x,2)+pow(y,2)+pow(z,2),0.5);}
	
	//------------COMPOUND ASSIGNMENT---------------
	vector& operator+=(const vector& rhs);
	vector& operator-=(const vector& rhs);
	vector& operator =(const vector& rhs);

	vector& operator*=(double rhs);	
	vector& operator/=(double rhs);	

	//-----------BINARY ARITHMETIC OPERATORS----------      
	const vector operator+(const vector& rhs) const ;
	const vector operator-(const vector& rhs) const ;
	
	const double operator*(const vector& rhs) const ; // dot product
	
	const vector operator*(double rhs) const ;
        const vector operator/(double rhs) const ;

	const bool operator==(const vector& rhs) const ;

};

const vector operator*(double lhs, const vector& rhs);

#endif
