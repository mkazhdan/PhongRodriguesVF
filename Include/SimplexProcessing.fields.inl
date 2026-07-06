/*
Copyright (c) 2026, Michael Kazhdan and Hongyi Liu
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

//////////////////////
// DerivativeTester //
//////////////////////
template< unsigned int K , HasDotProduct T >
template< HasSimplexFunctionAndFunctionDifferential< K , T > Field >
double DerivativeTester< K , T >::SquareError( Position< K > p , const Field & field , double delta )
{
	double error = 0;
	for( unsigned int k=0 ; k<K ; k++ )
	{
		// The offset positions along the k-th coordinate
		Point< double , K > _p[] = { p , p };
		_p[0][k] -= delta , _p[1][k] += delta;
		T d = ( field( _p[1] ) - field( _p[0] ) ) / ( 2. * delta ) - field.d( p )[k];
		error += DotProduct( d , d );
	}
	return error / K;
}

template< unsigned int K , HasDotProduct T >
template< HasSimplexFunctionAndFunctionDifferential< K , T > Field >
double DerivativeTester< K , T >::SquareError( const Field & field , unsigned int testCount , double delta )
{
	const Simplex< double , K , K > UnitRightSimplex = Simplex< double , K , K >::UnitRight();
	double error = 0;
	for( unsigned int c=0 ; c<testCount ; c++ ) error += SquareError( UnitRightSimplex.randomSample() , field , delta );
	return error / testCount;
}

////////////////////////
// NormalizationField //
////////////////////////
template< unsigned int K , HasDotProduct T , HasSimplexFunction< K , T > Field >
Differential< K , T > NormalizationField< K , T , Field >::d( Position< K > p ) const
{
	// N(x) = F(x) / < F(x) , F(x) >^0.5
	// dN(x) = dF(x) / < F(x) , F(x) >^0.5 - 0.5 * F(x) / < F(x) , F(x) >^1.5 * 2 * dF(x)
	//       = dF(x) / < F(x) , F(x) >^0.5 - F(x) / < F(x) , F(x) >^1.5 * dF(x)
	T f = _f(p);
	SimplexProcessing::Differential< K , T > df = _f.d(p);
	double l = sqrt( DotProduct( f , f ) );
	f /= l;
	for( unsigned int k=0 ; k<K ; k++ )
	{
		// Normalize
		df[k] /= l;
		// Project out the component in direction f
		df[k] -= DotProduct( df[k] , f ) * f;
	}
	return df;
}

///////////////////////
// LinearInterpolant //
///////////////////////
template< unsigned int K , typename T >
LinearInterpolant< K , T >::LinearInterpolant( void ){ for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = T(); }

template< unsigned int K , typename T >
LinearInterpolant< K , T >::LinearInterpolant( const T x[K+1] ){ for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = x[k]; }

template< unsigned int K , typename T >
T LinearInterpolant< K , T >::Value( Position< K > p , const T x[K+1] )
{
	T _x{};
	double p0 = 1.;
	for( unsigned int k=0 ; k<K ; k++ ) _x += x[k+1] * p[k] , p0 -= p[k];
	return _x + x[0] * p0;
}

template< unsigned int K , typename T >
Differential< K , T > LinearInterpolant< K , T >::DValue( Position< K > , const T x[K+1] )
{
	SimplexProcessing::Differential< K , T > d = {};
	for( unsigned int k=0 ; k<K ; k++ ) d[k] = x[k+1] - x[0];
	return d;
}

///////////////////////////////
// PhongRodriguesVectorField //
///////////////////////////////
template< unsigned int K , unsigned int N , bool Modulate >
typename PhongRodriguesVectorField< K , N , Modulate >::T PhongRodriguesVectorField< K , N , Modulate >::Value( Position< K > p , const T n[K+1] , const T x[K+1] )
{
	// Compute the normal at the inteprolated position
	T _n = LinearInterpolant< K , T >::Value( p , n );

	// Evaluated the Rodrigues rotation applied to the corner values
	T _x[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = SimplexProcessing::RodriguesRotation( n[k] , _n ) * x[k];

	if constexpr( Modulate ) return LinearInterpolant< K , T >::Value( p , _x );
	else
	{
		T sum = {};
		for( unsigned int k=0 ; k<=K ; k++ ) sum += _x[k];
		return sum;
	}
}

template< unsigned int K , unsigned int N , bool Modulate >
SimplexProcessing::Differential< K , typename PhongRodriguesVectorField< K , N , Modulate >::T > PhongRodriguesVectorField< K , N , Modulate >::DValue( Position< K > p , const T n[K+1] , const T x[K+1] )
{
	// F(p) = G( n(p) / |n(p)| )
	// dF = dG( n(p) / |n(p)| ) * d( n(p) / |n(p)| )
	// dF = dG( n(p) / |n(p)| ) * d( n(p) * < n(p) , n(p) >^{-1/2} )
	//    = dG( n(p) / |n(p)| ) * ( dn(p) * < n(p) , n(p) >^{-1/2} - 1/2 n(p) * < n(p) , n(p) >^{-3/2} * 2 * n(p) * dn(p) )
	//    = dG( n(p) / |n(p)| ) * ( dn(p) * < n(p) , n(p) >^{-1/2} - < n(p) , n(p) >^{-3/2} ( n(p) x n(p) ) * dn(p) )
	//    = dG( n(p) / |n(p)| ) * ( dn(p) - < n(p) , n(p) >^-1 ( n(p) x n(p) ) * dn(p) ) * < n(p) , n(p) >^{-1/2}
	//    = dG( n(p) / |n(p)| ) * ( dn(p) - ( n(p) x n(p) ) * dn(p) / |n(p)|^2 ) / | n(p) |

	// The interpolated normal
	T _n = LinearInterpolant< K , T >::Value( p , n );

	// R^k -> R^n
	// The differential of the interpolated normal
	SimplexProcessing::Differential< K , T > dn = LinearInterpolant< K , T >::DValue( p , n );

	T _x[K+1];
	SimplexProcessing::Differential< K , T > dX[K+1];
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		// Bring the tangent vector at the corner to the center point
		_x[k] = RodriguesRotation( n[k] , _n ) * x[k];

		// The differential of the rodrigues rotation with respect to the second argument, applied to the k-th tangent vector
		SquareMatrix< double , N > _dx;
		{
			SimplexProcessing::Differential< N , SquareMatrix< double , N > > dR = D2RodriguesRotation( n[k] , _n );
			for( unsigned int n=0 ; n<N ; n++ )
			{
				// The differential with respect to the n-th coordinate, applied to the k-th tangent vector
				Point< double , N > dR_x = dR[n] * x[k];
				for( unsigned int _n=0 ; _n<N ; _n++ ) _dx(n,_n) = dR_x[_n];
			}
		}
		for( unsigned int _k=0 ; _k<K ; _k++ ) dX[k][_k] = _dx * dn[_k];
	}

	if constexpr( Modulate ) return LinearInterpolant< K , T >::DValue( p , _x ) + LinearInterpolant< K , SimplexProcessing::Differential< K , T > >::Value( p , dX );
	else
	{
		SimplexProcessing::Differential< K , T > d;
		for( unsigned int k=0 ; k<=K ; k++ ) d += dX[k];
		return d;
	}
}

/////////////////////////////////////////////////////////
// PhongRodriguesIntrinsicToExtrinsicTangentXFormField //
/////////////////////////////////////////////////////////
template< unsigned int K >
PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Point< double , Dim > v[K+1] , const Point< double , Dim > n[K+1] )
	: _N(n)
{
	Point< double , Dim > tangents[K];
	for( unsigned int k=0 ; k<K ; k++ ) tangents[k] = v[k+1] - v[0];
	_normal = Point< double , Dim >::CrossProduct( tangents );
	for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _xForm(k,d) = tangents[k][d];
}

template< unsigned int K >
PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals )
	: PhongRodriguesIntrinsicToExtrinsicTangentXFormField( &vertices[0] , &normals[0] )
{}

template< unsigned int K >
typename PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::T PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::operator()( Position< K > p ) const
{
	return SimplexProcessing::RodriguesRotation( _normal , _N( p ) ) * _xForm;
}

template< unsigned int K >
SimplexProcessing::Differential< K , typename PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::T > PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::d( Position< K > p ) const
{
	SimplexProcessing::Differential< Dim , SquareMatrix< double , Dim > > dR = D2RodriguesRotation( _normal , _N(p) );
	SimplexProcessing::Differential< K , Point< double , Dim > > dN = _N.d(p);

	SimplexProcessing::Differential< K , Matrix< double , K , Dim > > d;
	for( unsigned int k=0 ; k<K ; k++ )
	{
		SquareMatrix< double , Dim > m;
		for( unsigned int d=0 ; d<Dim ; d++ ) m += dR[d] * dN[k][d];
		d[k] = m * _xForm;
	}

	return d;
}

/////////////////////////////////////////////////////////
// PhongRodriguesExtrinsicToIntrinsicTangentXFormField //
/////////////////////////////////////////////////////////
template< unsigned int K >
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Point< double , Dim > v[K+1] , const Point< double , Dim > n[K+1] )
	: _i2e( v , n ) , _gInv( MetricTensorFromEmbedding< K >( v ).inverse() )
{}

template< unsigned int K >
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals )
	: PhongRodriguesExtrinsicToIntrinsicTangentXFormField( &vertices[0] , &normals[0] )
{}

template< unsigned int K >
typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::operator()( Position< K > p ) const
{
	return _gInv * _i2e( p ).transpose();
}

template< unsigned int K >
SimplexProcessing::Differential< K , typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T > PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::d( Position< K > p ) const
{
	SimplexProcessing::Differential< K , typename PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::T > d = _i2e.d( p );
	SimplexProcessing::Differential< K , typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T > d_transpose;
	for( unsigned int k=0 ; k<K ; k++ ) d_transpose[k] = _gInv * d[k].transpose();
	return d_transpose;
}

////////////////////////////////
// ConnectionCoefficientField //
////////////////////////////////

template< unsigned int K >
ConnectionCoefficientField< K >::ConnectionCoefficientField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] ) 
	: _i2e( vertices , normals ) , _gInv( MetricTensorFromEmbedding< K >( vertices ).inverse() )
{}

template< unsigned int K >
ConnectionCoefficientField< K >::ConnectionCoefficientField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals ) 
	: ConnectionCoefficientField( &vertices[0] , &normals[0] )
{}

template< unsigned int K >
typename ConnectionCoefficientField< K >::T ConnectionCoefficientField< K >::operator()( Position< K > p ) const
{
	T C;

	Matrix< double , K , Dim > i2e = _i2e( p );
	SimplexProcessing::Differential< K , Matrix< double , K , Dim > > di2e = _i2e.d( p );

	// The coordinate acting as the vector field, and the coordinate along which we differentiate:
	for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ )
	{
		// The component of the derivative
		Point< double , K > dot;
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) dot[k] += i2e(k,d) * di2e[j](i,d);
		Point< double , K > coeff = _gInv * dot;
		for( unsigned int k=0 ; k<K ; k++ ) C(i,j,k) = coeff[k];
	}

	return C;
}

//////////////////////////
// IntrinsicVectorField //
//////////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename IntrinsicVectorField< K , N , VectorField >::T IntrinsicVectorField< K , N , VectorField >::operator()( Position< K > p ) const
{
	return _e2i(p) * _vf(p);
}

template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
SimplexProcessing::Differential< K , typename IntrinsicVectorField< K , N , VectorField >::T > IntrinsicVectorField< K , N , VectorField >::d( Position< K > p ) const
{
	SimplexProcessing::Differential< K , T > d;

	Point< double , N > v = _vf(p);
	SimplexProcessing::Differential< K , Point< double , N > > dv = _vf.d(p);

	Matrix< double , N , K > e2i = _e2i(p);
	SimplexProcessing::Differential< K , Matrix< double , N , K > > de2i = _e2i.d(p);

	for( unsigned int k=0 ; k<K ; k++ ) d[k] = de2i[k] * v + e2i * dv[k];

	return d;
}

////////////////////////////////
// SecondFundamentalFormField //
////////////////////////////////
template< unsigned int K , bool DifferentiateNormals >
SecondFundamentalFormField< K , DifferentiateNormals >::SecondFundamentalFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] ) 
	: _i2e( vertices , normals ) , _normals(normals)
{}

template< unsigned int K , bool DifferentiateNormals >
SecondFundamentalFormField< K , DifferentiateNormals >::SecondFundamentalFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals ) 
	: SecondFundamentalFormField( &vertices[0] , &normals[0] )
{}

template< unsigned int K , bool DifferentiateNormals >
typename SecondFundamentalFormField< K , DifferentiateNormals >::T SecondFundamentalFormField< K , DifferentiateNormals >::operator()( Position< K > p ) const
{
	if constexpr( DifferentiateNormals )
	{
		Differential< K , Point< double , Dim > > _dN = _normals.d(p);
		Matrix< double , K , Dim > dN;
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int n=0 ; n<K+1 ; n++ ) dN(k,n) = _dN[k][n];
		return _i2e(p).transpose() * dN;
	}
	else
	{
		T cov;
		Point< double , Dim > n = _normals(p);
		SimplexProcessing::Differential< K , Matrix< double , K , K+1 > > di2e = _i2e.d(p);

		for( unsigned int k=0 ; k<K ; k++ )
		{
			Point< double , K > e;
			e[k] = 1;
			for( unsigned int _k=0 ; _k<K ; _k++ ) cov(_k,k) = -Point< double , K+1 >::Dot( n , di2e[_k] * e );
		}
		return cov;
	}
}


//////////////////////////////
// CovariantDerivativeField //
//////////////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename CovariantDerivativeField< K , N , VectorField >::T
CovariantDerivativeField< K , N , VectorField >::operator()( Position< K > p ) const
{
	SquareMatrix< double , K > t;

	// The transformation from extrinsic to intrinsic tangents
	Matrix< double , N , K > e2i = _e2i(p);

	// The differentials of the extrinsic tangent vector fields
	Differential< K , Point< double , N > > dvf = _vf.d(p);

	for( unsigned int k=0 ; k<K ; k++ )
	{
		Point< double , K > _dx = e2i * dvf[k];
		for( unsigned int _k=0 ; _k<K ; _k++ ) t(k,_k) = _dx[_k];
	}

	return t;
}

/////////////////////////////////////////
// CovariantDirectionalDerivativeField //
/////////////////////////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunction< K , Point< double , N > > DirectionField , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename CovariantDirectionalDerivativeField< K , N , DirectionField , VectorField >::T
CovariantDirectionalDerivativeField< K , N , DirectionField , VectorField >::operator()( Position< K > p ) const
{
	// The transformation from extrinsic to intrinsic tangents
	Matrix< double , N , K > e2i = _e2i(p);

	// The evaluations of the vector fields, in the triangle tangent basis
	Point< double , K > dir = e2i( _dir(p) );

	// The differentials of the extrinsic tangent vector fields
	SimplexProcessing::Differential< K , Point< double , N > > dvf = _vf.d(p);

	// The covariant derivative
	return dvf( dir );
}

/////////////////////
// DivergenceField //
/////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename DivergenceField< K , N , VectorField >::T
DivergenceField< K , N , VectorField >::operator()( Position< K > p ) const
{
	return _dvf(p).trace();
}

//////////////////
// BracketField //
//////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField1 , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField2 >
typename BracketField< K , N , VectorField1 , VectorField2 >::T
BracketField< K , N , VectorField1 , VectorField2 >::operator()( Position< K > p ) const
{
	// The transformation from extrinsic to intrinsic tangents
	Matrix< double , N , K > e2i = _e2i(p);

	// The evaluations of the vector fields, in the triangle tangent basis
	Point< double , K > vf1 = e2i( _vf1(p) ) , vf2 = e2i( _vf2(p) );

	// The differentials of the extrinsic tangent vector fields
	Differential< K , Point< double , K+1 > > dvf1 = _vf1.d(p) , dvf2 = _vf2.d(p);

	// The bracket
	return dvf2( vf1 ) - dvf1( vf2 );
}