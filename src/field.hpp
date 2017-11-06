/*----------------------------------------------------------------------
 *  
 *   	FIPI - Fast Interface Particle Interactions
 *    
 *     	https://bottogroup.wordpress.com/
 *    
 *     	Chuan Gu, c.gu@qmul.ac.uk
 *    
 *     	Copyright (2016) Botto Research Group
 * 
 *	This software is distributed under the GNU General Public License.
 *  
 * ------------------------------------------------------------------------- */



#include"primitive.hpp"
#ifndef field_inc 
#define field_inc
//*********
class field
//*********
{
	private:
		fftw_complex** v; 
		std::string type;				
		int nx;
		int ny;
		int nz;
		int dimension; //1 for scalar 3 for vector 9 for tensor
		bool in_fspace;
	public:

		//-------------CONSTRUCTORS AND DESTRUCTORS-----------
		field(std::string type_id);		
		field(std::string type_id, std::string name_id);
		~field();	
		
		void malloc(std::string space);
		void reset(std::string space);
		void set_value (int d, int it, double value);	
		void add_value (int d, int it, double value);	
		//!		
		const field norm() const;
		double interpolate (const vector& position, int d) const;
		double interpolate2(const vector& position, int d) const;

		//------------VECTOR CALCULUS-----------
		const field grad(const field& kf) const;	
		const field div (const field& kf) const;
		const field dot  (const field& b) const;
		const field power(double b) const;

		//------------FFT----------------------	
		void fft (const fftw_plan& p_forward);
		void ifft(const fftw_plan& p_backward);
		void free();
	
		//------------OUTPUT--------------------	
		void output_to_file(std::string file_name, std::string variable_name);

		//------------COMPOUND ASSIGNMENT---------------
		field& operator+=(const field& rhs);
		field& operator-=(const field& rhs);
		field& operator*=(const field& rhs);
		field& operator/=(const field& rhs);

		field& operator+=(double rhs);
		field& operator-=(double rhs);
		field& operator*=(double rhs);
		field& operator/=(double rhs);
		
		field& operator=(const field& b);
	
		//-----------BINARY ARITHMETIC OPERATORS----------	
		const field operator+(const field& rhs) const ;
		const field operator-(const field& rhs) const ;
		const field operator*(const field& rhs) const ;	
		const field operator/(const field& rhs) const ;
		
		const field operator+(double rhs) const ;
		const field operator-(double rhs) const ;
		const field operator*(double rhs) const ;
		const field operator/(double rhs) const ;

		void print(int it);		
};

const field operator+(double b, const field& a);

const field operator*(double b, const field& a);

const field operator-(double b, const field& a);
#endif
