// _________________________________________________________________________
//
//  TEMPLATE CLASS IMPLEMENTATION OF TRUNCATED POWER SERIES AND ITS ALGEBRA
// _________________________________________________________________________
//
// TPSA
// ====
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		August 31st, 2006 @ LNLS
//
// Notes:
//
// 01. This is an implementation for a pre-defined-at-compile-time number of variables
//     and of polynomial order.
// 02. Multiplication algorithm makes use of one OSIP (one-step index pointer) in a
//     forward scheme (as thought by A. Dragt, according to Y. Yan) in order to achieve
//     optimized efficiency.
// 03. I have to better understand template friend functions and operators. This way I
//     should be able to convert some non-member functions to friends and improve
//     their algorithms by having access to protected class members. Try to ask folks
//     from the C++ community?
// 04. Is it worth trying to implement syntactic sugars with expression templates in order
//     to minimize the problem of temporaries? It is not clear for me that, for the TPS class,
//     something is to be gained in terms of efficiency...
// 05. Is it worth trying to implement multiplication with FFT? How to map multivariate
//     polynomials into (can it be onto?) univariate polynomials of higher order? The prospects
//     of using FFT are interesting since the Beam Dynamics community does not seem to have
//     realized that the method could speed up calculation of taylor maps...



#ifndef TPSA_H
#define TPSA_H

#include <complex>
#include <cstring>
#include <ostream>


// Expression Templates: IMPLEMENTATION OF BINOMIALS COEFFICIENTS AND RELATED RELEVANT EXPRESSIONS
// -----------------------------------------------------------------------------------------------

template <int N, int K>        struct et_binomial      { enum { val = et_binomial<N-1,K-1>::val + et_binomial<N-1,K>::val }; };
template <int K>               struct et_binomial<0,K> { enum { val = 0 }; };
template <>                    struct et_binomial<0,0> { enum { val = 1 }; };
template <int V, int N, int n> struct et_osip          { enum { val = (V * et_binomial<V+n,n>::val * et_binomial<V+N-n,N-n>::val)/(V+n) + et_osip<V,N,n-1>::val }; };
template <int V, int N>        struct et_osip<V,N,0>   { enum { val = et_binomial<V,0>::val * et_binomial<V+N,N>::val }; };


// Forward Declarations
// --------------------

template <unsigned int  , unsigned int  , typename     > class Tpsa;
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> abs   (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> sqrt  (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> log   (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> atan  (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> atanh (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> cos   (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> sin   (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> tan   (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> cosh  (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> sinh  (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> asin  (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> acos  (const Tpsa<V,N,TYPE>&);
template <unsigned int V, unsigned int N, typename TYPE> Tpsa<V,N,TYPE> D     (const Tpsa<V,N,TYPE>&, const unsigned int);


// Definition: CLASS TPS
// ---------------------

template <unsigned int V = 0, unsigned int N = 0, typename TYPE = double>
class Tpsa {

	friend Tpsa<V,N,TYPE> abs<>   (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> sqrt<>  (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> log<>   (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> atan<>  (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> atanh<> (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> cos<>   (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> sin<>   (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> tan<>   (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> cosh<>  (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> sinh<>  (const Tpsa<V,N,TYPE>&);
	template <typename T> friend Tpsa<V,N,T> asin (const Tpsa<V,N,T>&);
	template <typename T> friend Tpsa<V,N,T> acos (const Tpsa<V,N,T>&);
	//friend Tpsa<V,N,TYPE> asin<>  (const Tpsa<V,N,TYPE>&);
	friend Tpsa<V,N,TYPE> D<>     (const Tpsa<V,N,TYPE>&, const unsigned int);

public:

	// constructors
	Tpsa(const TYPE& a_ = 0, const unsigned int v_ = V);
	Tpsa(const Tpsa& a_);

	// algebra of class elements with support field
	Tpsa  operator *  (const TYPE& o_) const;
	Tpsa  operator /  (const TYPE& o_) const;
	Tpsa& operator += (const TYPE& o_);
	Tpsa& operator -= (const TYPE& o_);
	Tpsa& operator *= (const TYPE& o_);
	Tpsa& operator /= (const TYPE& o_);

	// algebra of class elements
	Tpsa  operator -  () const;
	Tpsa  operator +  (const Tpsa& o_) const;
	Tpsa  operator -  (const Tpsa& o_) const;
	Tpsa  operator *  (const Tpsa& o_) const;
	Tpsa  inverse     ()              const;
	Tpsa  operator /  (const Tpsa& o_) const;
	Tpsa& operator += (const Tpsa& o_);
	Tpsa& operator -= (const Tpsa& o_);
	Tpsa& operator *= (const Tpsa& o_);
	Tpsa& operator /= (const Tpsa& o_);

	// boolean operators
	bool operator == (const TYPE& o_) const;
	bool operator != (const TYPE& o_) const;
	bool operator <  (const TYPE& o_) const;
	bool operator <= (const TYPE& o_) const;
	bool operator >  (const TYPE& o_) const;
	bool operator >= (const TYPE& o_) const;

	bool operator >  (const Tpsa& o_) const;
	bool operator >= (const Tpsa& o_) const;
	bool operator == (const Tpsa& o_) const;
	bool operator != (const Tpsa& o_) const;

	explicit operator int()    const { return int(c[0]); }
	explicit operator double() const { return double(c[0]); }

	// auxiliary public functions
	static unsigned int  get_n()         { return N; }
	static unsigned int  get_v()         { return V; }
	static unsigned int  get_size()      { return et_binomial<N+V,V>::val; }
	static unsigned int  get_osip_size() { return et_osip<V,N,N>::val; }
	static unsigned int  get_index (const unsigned int* power_);
	static void          get_power (const unsigned int idx, unsigned int* power);
	const TYPE&          get_c(unsigned int index) const { return c[index]; }
	TYPE&                set_c(unsigned int index)  { return c[index]; }

//private:
public:

	TYPE                c[et_binomial<N+V,V>::val];
	static bool         init;
	static unsigned int binomials[((N+V+1)*(N+V+2))>>1];
	static void         initialization();
	static unsigned int osip[et_osip<V,N,N>::val];
	static unsigned int powers[et_osip<V,N,N>::val][V];

	static unsigned int& C(unsigned int v, unsigned int n) { return binomials[(((v+n-1)*(n+v))>>1) + n]; }
	static unsigned int  first_at_order(unsigned int order) { return (order==0) ? 0 : C(V,order-1); }
	static unsigned int  last_at_order (unsigned int order) { return (order==0) ? 1 : C(V,order); }

};



// Static Members
// --------------

template <unsigned int V, unsigned int N, typename TYPE> bool Tpsa<V,N,TYPE>::init = true;
template <unsigned int V, unsigned int N, typename TYPE> unsigned int Tpsa<V,N,TYPE>::binomials[((N+V+1)*(N+V+2))>>1];
template <unsigned int V, unsigned int N, typename TYPE> unsigned int Tpsa<V,N,TYPE>::osip[et_osip<V,N,N>::val];
template <unsigned int V, unsigned int N, typename TYPE> unsigned int Tpsa<V,N,TYPE>::powers[et_osip<V,N,N>::val][V];


// Implementations: CONSTRUCTORS
// -----------------------------

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>::Tpsa(const TYPE& a_, const unsigned int v_) {
	if (init) initialization();
	memset(this->c, 0, sizeof(TYPE)*get_size());
	//for(unsigned int i=1; i<get_size(); i++) c[i] = 0;
	c[0] = a_;
	if (v_<V) c[v_+1] = 1;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>::Tpsa(const Tpsa<V,N,TYPE>& a_)  {
	for(unsigned int i=0; i<get_size(); i++) c[i] = a_.c[i];
}



// Implementations: ALGEBRA OF CLASS ELEMENTS WITH FIELD
// -----------------------------------------------------

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator * (const TYPE& o_) const {
	Tpsa<V,N,TYPE> r(*this);
	for(unsigned int i=0; i<get_size(); i++) r.c[i] *= o_;
	return r;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator / (const TYPE& o_) const {
	Tpsa<V,N,TYPE> r(*this);
	for(unsigned int i=0; i<get_size(); i++) r.c[i] /= o_;
	return r;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator += (const TYPE& o_) {
	c[0] += o_;
	return *this;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator -= (const TYPE& o_) {
	c[0] -= o_;
	return *this;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator *= (const TYPE& o_) {
	for(unsigned int i=0; i<get_size(); i++) c[i] *= o_;
	return *this;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator /= (const TYPE& o_) {
	for(unsigned int i=0; i<get_size(); i++) c[i] /= o_;
	return *this;
};


// Implementations: ALGEBRA OF CLASS ELEMENTS
// ------------------------------------------

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator - () const {
	Tpsa<V,N,TYPE> r(*this);
	for(unsigned int i=0; i<get_size(); i++) r.c[i] = -r.c[i];
	return r;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator + (const Tpsa<V,N,TYPE>& o_) const {
	Tpsa<V,N,TYPE> r(*this);
	for(unsigned int i=0; i<get_size(); i++) r.c[i] += o_.c[i];
	return r;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator - (const Tpsa<V,N,TYPE>& o_) const {
	Tpsa<V,N,TYPE> r(*this);
	for(unsigned int i=0; i<get_size(); i++) r.c[i] -= o_.c[i];
	return r;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator * (const Tpsa<V,N,TYPE>& o_) const {
	Tpsa<V,N,TYPE> r;
	unsigned int counter = 0;
	for(unsigned int n1 = 0; n1 <= N; n1++) {
		for(unsigned int i1 = first_at_order(n1); i1 < last_at_order(n1); i1++) {
			for(unsigned int n2 = 0; n1+n2 <= N; n2++) {
				for(unsigned int i2 = first_at_order(n2); i2 < last_at_order(n2); i2++) {
					r.c[osip[counter]] += this->c[i1] * o_.c[i2];
					counter++;
				}
			}
		}
	}
	return r;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::inverse() const {
	Tpsa<V,N,TYPE> r;
	TYPE          a = this->c[0];
	Tpsa<V,N,TYPE> x(*this); x.c[0] = 0; x /= a;
	Tpsa<V,N,TYPE> p(1);
	for(unsigned int i=0; i<=N; i++) {
		r += p * (i&1?-1:1);
		p *= x;
	}
	return r / a;
}

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> Tpsa<V,N,TYPE>::operator / (const Tpsa<V,N,TYPE>& o_) const {
	return *this * o_.inverse();
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator += (const Tpsa<V,N,TYPE>& o_) {
	for(unsigned int i=0; i<get_size(); i++) c[i] += o_.c[i];
	return *this;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator -= (const Tpsa<V,N,TYPE>& o_) {
	for(unsigned int i=0; i<get_size(); i++) c[i] -= o_.c[i];
	return *this;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator *= (const Tpsa<V,N,TYPE>& o_) {
	*this = *this * o_;
	return *this;
};

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE>& Tpsa<V,N,TYPE>::operator /= (const Tpsa<V,N,TYPE>& o_) {
	*this = *this * o_.inverse();
	return *this;
};


// Implementation: BOOLEAN OPERATORS
// ---------------------------------

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator == (const TYPE& o_) const {
	return c[0] == o_;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator != (const TYPE& o_) const {
	return c[0] != o_;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator < (const TYPE& o_) const {
	for(unsigned int i=0; i<get_size(); i++) if (c[i] - o_ != 0) return c[i] - o_ < 0;
	return false;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator <= (const TYPE& o_) const {
	for(unsigned int i=0; i<get_size(); i++) if (c[i] - o_ != 0) return c[i] - o_ <= 0;
	return false;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator > (const TYPE& o_) const {
	for(unsigned int i=0; i<get_size(); i++) if (c[i] - o_ != (TYPE) 0) return c[i] - o_ > 0;
	return false;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator >= (const TYPE& o_) const {
	for(unsigned int i=0; i<get_size(); i++) if (c[i] - o_ != (TYPE) 0) return c[i] - o_ >= 0;
	return false;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator >= (const Tpsa<V,N,TYPE>& o_) const {
	return (*this - o_) >= (TYPE) 0;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator == (const Tpsa<V,N,TYPE>& o_) const {
	return (*this - o_) == (TYPE) 0;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator != (const Tpsa<V,N,TYPE>& o_) const {
	return (*this - o_) != (TYPE) 0;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool Tpsa<V,N,TYPE>::operator > (const Tpsa<V,N,TYPE>& o_) const {
	return (*this - o_) > (TYPE) 0;
}


// Implementation: AUXILIARY MEMBER FUNCTIONS
// ------------------------------------------

template <unsigned int V, unsigned int N, typename TYPE>
void Tpsa<V,N,TYPE>::initialization() {
	// sets binomial coefficients
	for(unsigned int s=0; s<=N+V; s++) {
		C(s,0) = C(0,s) = 1;
		for(unsigned int n=1; n < s; n++) C(s-n,n) = C((s-1)-n,n) + C((s-1)-(n-1),(n-1));
	}
	// sets one-step index pointers
	unsigned int counter = 0;
	unsigned int power1[V], power2[V], power[V];
	for(unsigned int n1 = 0; n1 <= N; n1++) {
		for(unsigned int i1 = first_at_order(n1); i1 < last_at_order(n1); i1++) {
			get_power(i1,power1);
			for(unsigned int n2 = 0; n1+n2 <= N; n2++) {
				for(unsigned int i2 = first_at_order(n2); i2 < last_at_order(n2); i2++) {
					get_power(i2,power2);
					for(unsigned int i=0; i<V; i++) {
						power[i] = power1[i] + power2[i];
						powers[counter][i] = power[i];
					}
					osip[counter] = get_index(power);
					counter++;
				}
			}
		}
	}
	init = false;
}

template <unsigned int V, unsigned int N, typename TYPE>
unsigned int Tpsa<V,N,TYPE>::get_index(const unsigned int* power_) {
	unsigned int s_i = 0;
	unsigned int idx = 0;
	for(unsigned int i=0; i<V-1; i++) {
		s_i += power_[V-1-i];
		idx += (s_i * C(i,s_i)) / (i+1);
	}
	s_i += power_[0];
	idx = C(V,s_i) - idx;
	return idx - 1;
}

template <unsigned int V, unsigned int N, typename TYPE>
void Tpsa<V,N,TYPE>::get_power(unsigned int idx_, unsigned int* power_) {
	unsigned int idx = idx_ + 1;
	for(power_[0]=0; C(V,power_[0]) < idx; power_[0]++);
	idx = C(V,power_[0]) - idx;
	for(unsigned int i=1; i<V; i++) {
		for(power_[i]=power_[i-1]; (power_[i]*C(V-i-1,power_[i]))/(V-i) > idx; power_[i]--);
		idx -= (power_[i]*C(V-i-1,power_[i]))/(V-i);
		power_[i-1] -= power_[i];
	}
}


// Implementation: NON-MEMBER FUNCTIONS AND OPERATORS WITH ARGUMENTS OF CLASS TYPE
// -------------------------------------------------------------------------------
// Note: these functions are generating Stack Overflow at run-time for large N,V
// (N=10,V=10 for example) at a P4 with 1Gb of RAM

template <typename T, unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> operator + (const T& o1, const Tpsa<V,N,TYPE>& o2) { return o2 + o1; }

template <typename T, unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> operator * (const T& o1, const Tpsa<V,N,TYPE>& o2) { return o2 * o1; }

template <typename T, unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> operator - (const T& o1, const Tpsa<V,N,TYPE>& o2) { return (-o2) + o1; }

template <typename T, unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> operator / (const T& o1, const Tpsa<V,N,TYPE>& o2) { return o2.inverse() * o1; }


template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> abs(const Tpsa<V,N,TYPE>& a_) {
	if (a_ >= 0) return a_; else return -a_;
}

//template <unsigned int V, unsigned int N, typename TYPE>
//Tpsa<V,N,TYPE> fabs(const Tpsa<V,N,TYPE>& a_) { return abs(a_); }

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> sqrt(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> r;
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0; x /= a_.c[0];
	Tpsa<V,N,TYPE> p(1);
	TYPE          f = 1;
	for(unsigned int i=0; i<=N; i++) {
		r += p * f;
		f *= (0.5 - i)/(i+1);
		p *= x;
	}
	r *= sqrt(a_.c[0]);
	return r;
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> log(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> r;
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0; x /= a_.c[0];;
	Tpsa<V,N,TYPE> p(x);
	for(int i=1; i<=N; i++) {
		r += p / ((i&1?1.0:-1.0)*i);
		p *= x;
	}
	r += log(a_.c[0]);
	return r;
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> cos(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> rc(1), rs;
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0;
	Tpsa<V,N,TYPE> p(x);
	unsigned int fac = 1;
	for(int i=1; i<=N; i++) {
		if (i&1) { rs += ((i&2)?-1:1) * p / fac; } else { rc += ((i&2)?-1:1) * p / fac; }
		p   *= x;
		fac *= i+1;
	}
	return cos(a_.c[0]) * rc - sin(a_.c[0]) * rs;
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> sin(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> rc(1), rs;
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0;
	Tpsa<V,N,TYPE> p(x);
	unsigned int fac = 1;
	for(int i=1; i<=N; i++) {
		if (i&1) { rs += ((i&2)?-1:1) * p / fac; } else { rc += ((i&2)?-1:1) * p / fac; }
		p   *= x;
		fac *= i+1;
	}
	return sin(a_.c[0]) * rc + cos(a_.c[0]) * rs;
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> tan(const Tpsa<V,N,TYPE>& a_) {
	return sin(a_)/cos(a_);
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> cosh(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> rc(1), rs;
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0;
	Tpsa<V,N,TYPE> p(x);
	unsigned int fac = 1;
	for(int i=1; i<=N; i++) {
		if (i&1) { rs += p / fac; } else { rc += p / fac; }
		p   *= x;
		fac *= i+1;
	}
	return cosh(a_.c[0]) * rc + sinh(a_.c[0]) * rs;
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> sinh(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> rc(1), rs;
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0;
	Tpsa<V,N,TYPE> p(x);
	unsigned int fac = 1;
	for(int i=1; i<=N; i++) {
		if (i&1) { rs += p / fac; } else { rc += p / fac; }
		p   *= x;
		fac *= i+1;
	}
	return sinh(a_.c[0]) * rc + cosh(a_.c[0]) * rs;
}


#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> atan(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> r(atan(a_.c[0]));
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0;
	Tpsa<V,N,TYPE> p(x);
	TYPE          q  = sqrt(1 + a_.c[0]*a_.c[0]);
	TYPE          pq = q;
	TYPE          c1 = 1/q;
	TYPE          s1 = a_.c[0] * c1;
	TYPE          ci = c1;
	TYPE          si = s1;
	// ArcTan function implementation based on expression D[1/(1+x^2)] = (1/2) D[1/(1+xi) + 1/(1-xi)]
	for(int i=1; i<=N; i++) {
		if (i&1) {
			// odd terms
			r += p * (i&2?-1:1)*ci/(i*pq);
		} else {
			// even terms
			r += p * (i&2?-1:1)*si/(i*pq);
		}
		TYPE ts = s1 * ci + c1 * si;
		ci = c1 * ci - s1 * si;
		si = ts;
		pq *= q;
		p  *= x;
	}
	return r;
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> asin(const Tpsa<V,N,TYPE>& a_) {

	std::complex<TYPE> I(0,1);
	Tpsa<V,N,TYPE> t1 = sqrt(1-a_*a_);
	Tpsa<V,N,std::complex<TYPE> > t2;
	for(unsigned int i=0; i<a_.get_size(); i++) t2.c[i] = I * a_.c[i] + t1.c[i];
	Tpsa<V,N,std::complex<TYPE> > t3 = log(t2);
	Tpsa<V,N,TYPE> r;
	for(unsigned int i=0; i<a_.get_size(); i++) r.c[i] = t3.c[i].imag();
	return r;

}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> acos(const Tpsa<V,N,TYPE>& a_) {
	return M_PI/2 - asin(a_);
}

#include <math.h>
template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> atanh(const Tpsa<V,N,TYPE>& a_) {
	Tpsa<V,N,TYPE> r(atan(a_.c[0]));
	Tpsa<V,N,TYPE> x(a_); x.c[0] = 0;
	Tpsa<V,N,TYPE> p(x);
	TYPE          q  = sqrt(1 + a_.c[0]*a_.c[0]);
	TYPE          pq = q;
	TYPE          c1 = 1/q;
	TYPE          s1 = a_.c[0] * c1;
	TYPE          ci = c1;
	TYPE          si = s1;
	// ArcTan function implementation based on expression D[1/(1+x^2)] = (1/2) D[1/(1+xi) + 1/(1-xi)]
	for(int i=1; i<=N; i++) {
		if (i&1) {
			// odd terms
			r += p * (i&2?-1:1)*ci/(i*pq);
		} else {
			// even terms
			r += p * (i&2?-1:1)*si/(i*pq);
		}
		TYPE ts = s1 * ci + c1 * si;
		ci = c1 * ci - s1 * si;
		si = ts;
		pq *= q;
		p  *= x;
	}
	return r;
}

template <unsigned int V, unsigned int N, typename TYPE>
Tpsa<V,N,TYPE> D(const Tpsa<V,N,TYPE>& a_, const unsigned int v_) {
	if (v_>=V) return 0;
	Tpsa<V,N,TYPE> r;
	unsigned int power[V];
	for(unsigned int i=0; i<a_.get_size(); i++) {
		a_.get_power(i,power);
		unsigned int f = power[V-1-v_];
		if (f!=0) {
			power[V-1-v_]--;
			unsigned int idx = a_.get_index(power);
			r.c[idx] += a_.c[i] * (TYPE) f;
		}
	}
	return r;
}

#include <iostream>
template <int V, int N, typename TYPE>
std::ostream& operator << (std::ostream& out, const Tpsa<V,N,TYPE>& o) {
	//out << Tpsa<V,N,TYPE>::get_v() << " " << Tpsa<V,N,TYPE>::get_n() << " " << Tpsa<V,N,TYPE>::get_size() << std::endl;
	unsigned int p[V];
	for(unsigned int i=0; i<Tpsa<V,N,TYPE>::get_size(); i++) {
		out << i << "  ";
		Tpsa<V,N,TYPE>::get_power(i,p);
		for(unsigned int v=0; v<Tpsa<V,N,TYPE>::get_v(); v++) out << p[v] << " ";
		out << "  " << o.get_c(i); out << std::endl;
	}
	return out;
}

template <unsigned int V, unsigned int N, typename TYPE>
bool isfinite(const Tpsa<V,N,TYPE>& a_) {
	return std::isfinite(a_.c[0]);
}




#endif
