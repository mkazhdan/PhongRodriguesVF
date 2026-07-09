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

///////////////////////
// LinearInterpolant //
///////////////////////
template< unsigned int K , typename T >
LinearInterpolant< K , T >::LinearInterpolant( void ){ for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = T(); }

template< unsigned int K , typename T >
LinearInterpolant< K , T >::LinearInterpolant( const T x[K+1] ){ for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = x[k]; }

template< unsigned int K , typename T >
T LinearInterpolant< K , T >::Value( const T x[K+1] , Position< K > p )
{
	T _x{};
	double p0 = 1.;
	for( unsigned int k=0 ; k<K ; k++ ) _x += x[k+1] * p[k] , p0 -= p[k];
	return _x + x[0] * p0;
}

template< unsigned int K , typename T >
Differential< K , T > LinearInterpolant< K , T >::DValue( const T x[K+1] , Position< K > )
{
	SimplexProcessing::Differential< K , T > d = {};
	for( unsigned int k=0 ; k<K ; k++ ) d[k] = x[k+1] - x[0];
	return d;
}

///////////////////////////////
// PhongRodriguesVectorField //
///////////////////////////////
template< unsigned int K , unsigned int N , bool Modulate >
typename PhongRodriguesVectorField< K , N , Modulate >::T PhongRodriguesVectorField< K , N , Modulate >::Value( const T n[K+1] , const T x[K+1] , Position< K > p )
{
	// Compute the normal at the inteprolated position
	T _n = LinearInterpolant< K , T >::Value( n , p );

	// Evaluated the Rodrigues rotation applied to the corner values
	T _x[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = SimplexProcessing::RodriguesRotation( n[k] , _n ) * x[k];

	if constexpr( Modulate ) return LinearInterpolant< K , T >::Value( _x , p );
	else
	{
		T sum = {};
		for( unsigned int k=0 ; k<=K ; k++ ) sum += _x[k];
		return sum;
	}
}

template< unsigned int K , unsigned int N , bool Modulate >
SimplexProcessing::Differential< K , typename PhongRodriguesVectorField< K , N , Modulate >::T > PhongRodriguesVectorField< K , N , Modulate >::DValue( const T n[K+1] , const T x[K+1] , Position< K > p )
{
	// F(p) = G( n(p) / |n(p)| )
	// dF = dG( n(p) / |n(p)| ) * d( n(p) / |n(p)| )
	// dF = dG( n(p) / |n(p)| ) * d( n(p) * < n(p) , n(p) >^{-1/2} )
	//    = dG( n(p) / |n(p)| ) * ( dn(p) * < n(p) , n(p) >^{-1/2} - 1/2 n(p) * < n(p) , n(p) >^{-3/2} * 2 * n(p) * dn(p) )
	//    = dG( n(p) / |n(p)| ) * ( dn(p) * < n(p) , n(p) >^{-1/2} - < n(p) , n(p) >^{-3/2} ( n(p) x n(p) ) * dn(p) )
	//    = dG( n(p) / |n(p)| ) * ( dn(p) - < n(p) , n(p) >^-1 ( n(p) x n(p) ) * dn(p) ) * < n(p) , n(p) >^{-1/2}
	//    = dG( n(p) / |n(p)| ) * ( dn(p) - ( n(p) x n(p) ) * dn(p) / |n(p)|^2 ) / | n(p) |

	// The interpolated normal
	T _n = LinearInterpolant< K , T >::Value( n , p );

	// R^k -> R^n
	// The differential of the interpolated normal
	SimplexProcessing::Differential< K , T > dn = LinearInterpolant< K , T >::DValue( n , p );

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

	if constexpr( Modulate ) return LinearInterpolant< K , T >::DValue( _x , p ) + LinearInterpolant< K , SimplexProcessing::Differential< K , T > >::Value( dX , p );
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
{
	Point< double , Dim > tangents[K];
	for( unsigned int k=0 ; k<=K ; k++ ) _normals[k] = n[k];
	for( unsigned int k=0 ; k<K ; k++ ) tangents[k] = v[k+1] - v[0];
	_normal = Point< double , Dim >::CrossProduct( tangents );
	for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _xForm(k,d) = tangents[k][d];
}

template< unsigned int K >
PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals )
	: PhongRodriguesIntrinsicToExtrinsicTangentXFormField( &vertices[0] , &normals[0] )
{}

template< unsigned int K >
typename PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::T
PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::Value
(
	const Point< double , Dim > & normal ,
	const Point< double , Dim > normals[K+1] ,
	const Matrix< double , K , Dim > & xForm ,
	Position< K > p
)
{
	return SimplexProcessing::RodriguesRotation( normal , LinearInterpolant< K , Point< double , Dim > >::Value( normals , p ) ) * xForm;
}

template< unsigned int K >
SimplexProcessing::Differential< K , typename PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::T >
PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::DValue
(
	const Point< double , Dim > & normal ,
	const Point< double , Dim > normals[K+1] ,
	const Matrix< double , K , Dim > & xForm ,
	Position< K > p
)
{
	SimplexProcessing::Differential< Dim , SquareMatrix< double , Dim > > dR = D2RodriguesRotation( normal , LinearInterpolant< K , Point< double , Dim > >::Value( normals , p ) );
	SimplexProcessing::Differential< K , Point< double , Dim > > dN = LinearInterpolant< K , Point< double , Dim > >::DValue( normals , p );

	SimplexProcessing::Differential< K , Matrix< double , K , Dim > > d;
	for( unsigned int k=0 ; k<K ; k++ )
	{
		SquareMatrix< double , Dim > m;
		for( unsigned int d=0 ; d<Dim ; d++ ) m += dR[d] * dN[k][d];
		d[k] = m * xForm;
	}

	return d;
}

/////////////////////////////////////////////////////////
// PhongRodriguesExtrinsicToIntrinsicTangentXFormField //
/////////////////////////////////////////////////////////
template< unsigned int K >
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Point< double , Dim > v[K+1] , const Point< double , Dim > n[K+1] )
	: PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >( v , n ) , _gInv( MetricTensorFromEmbedding< K >( v ).inverse() )
{}

template< unsigned int K >
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals )
	: PhongRodriguesExtrinsicToIntrinsicTangentXFormField( &vertices[0] , &normals[0] )
{}

template< unsigned int K >
typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Matrix< double , K , Dim > & i2e
)
{
	return gInv * i2e.transpose();
}

template< unsigned int K >
SimplexProcessing::Differential< K , typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T >
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::DValue
(
	const SquareMatrix< double , K > & gInv ,
	const SimplexProcessing::Differential< K , Matrix< double , K , Dim > > & di2e
)
{
	SimplexProcessing::Differential< K , typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T > d;
	for( unsigned int k=0 ; k<K ; k++ ) d[k] = gInv * di2e[k].transpose();
	return d;
}

template< unsigned int K >
typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , Dim > & normal ,
	const Point< double , Dim > normals[K+1] ,
	const Matrix< double , K , Dim > & xForm ,
	Position< K > p
)
{
	return Value( gInv , PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::Value( normal , normals , xForm , p ) );
}

template< unsigned int K >
SimplexProcessing::Differential< K , typename PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T >
PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::DValue
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , Dim > & normal ,
	const Point< double , Dim > normals[K+1] ,
	const Matrix< double , K , Dim > & xForm ,
	Position< K > p
)
{
	return DValue( gInv , PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::DValue( normal , normals , xForm , p ) );
}

////////////////////////////////
// ConnectionCoefficientField //
////////////////////////////////
template< unsigned int K , bool Symmetrize >
ConnectionCoefficientField< K , Symmetrize >::ConnectionCoefficientField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] ) 
	: PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >( vertices , normals ) , _gInv( MetricTensorFromEmbedding< K >( vertices ).inverse() )
{}

template< unsigned int K , bool Symmetrize >
ConnectionCoefficientField< K , Symmetrize >::ConnectionCoefficientField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals ) 
	: ConnectionCoefficientField( &vertices[0] , &normals[0] )
{}

template< unsigned int K , bool Symmetrize >
typename ConnectionCoefficientField< K , Symmetrize >::T
ConnectionCoefficientField< K , Symmetrize >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Matrix< double , K , Dim > & i2e ,
	const Differential< K , Matrix< double , K , Dim > > & di2e
)
{
	T C;

	// The coordinate acting as the vector field, and the coordinate along which we differentiate:
	for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ )
	{
		// The component of the derivative
		Point< double , K > dot;
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) dot[k] += i2e(k,d) * di2e[j](i,d);
		Point< double , K > coeff = gInv * dot;
		for( unsigned int k=0 ; k<K ; k++ ) C(i,j,k) = coeff[k];
	}

	if constexpr( Symmetrize )
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int l=0 ; l<k ; l++ ) for( unsigned int m=0 ; m<K ; m++ ) C(k,l,m) = C(l,k,m) = ( C(k,l,m) + C(l,k,m) ) / 2.;

	return C;
}

template< unsigned int K , bool Symmetrize >
typename ConnectionCoefficientField< K , Symmetrize >::T
ConnectionCoefficientField< K , Symmetrize >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , Dim > & normal ,
	const Point< double , Dim > normals[K+1] ,
	const Matrix< double , K , Dim > & xForm ,
	Position< K > p
)
{
	return Value( gInv , PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::Value( normal , normals , xForm , p ) , PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::DValue( normal , normals , xForm , p ) );
}

//////////////////////////
// IntrinsicVectorField //
//////////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename IntrinsicVectorField< K , N , VectorField >::T
IntrinsicVectorField< K , N , VectorField >::Value
(
	const Matrix< double , N , K > & e2i ,
	const Point< double , N > & v
)
{
	return e2i * v;
}

template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
SimplexProcessing::Differential< K , typename IntrinsicVectorField< K , N , VectorField >::T >
IntrinsicVectorField< K , N , VectorField >::DValue
(
	const Matrix< double , N , K > & e2i ,
	const SimplexProcessing::Differential< K , Matrix< double , N , K > > & de2i ,
	const Point< double , N > & v ,
	const SimplexProcessing::Differential< K , Point< double , N > > & dv
)
{
	SimplexProcessing::Differential< K , T > d;
	for( unsigned int k=0 ; k<K ; k++ ) d[k] = de2i[k] * v + e2i * dv[k];
	return d;
}

template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename IntrinsicVectorField< K , N , VectorField >::T
IntrinsicVectorField< K , N , VectorField >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const VectorField & VF ,
	Position< K > p
)
{
	return Value( PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value( gInv , normal , normals , xForm , p ) , VF(p) );
}

template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
SimplexProcessing::Differential< K , typename IntrinsicVectorField< K , N , VectorField >::T >
IntrinsicVectorField< K , N , VectorField >::DValue
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const VectorField & VF ,
	Position< K > p
)
{
	return DValue( PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value( gInv , normal , normals , xForm , p ) , PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::DValue( gInv , normal , normals , xForm , p ) , VF(p) , VF.d(p) );
}

////////////////////////////////
// SecondFundamentalFormField //
////////////////////////////////
template< unsigned int K , bool DifferentiateNormals >
SecondFundamentalFormField< K , DifferentiateNormals >::SecondFundamentalFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] ) 
	: PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >( vertices , normals )
{}

template< unsigned int K , bool DifferentiateNormals >
SecondFundamentalFormField< K , DifferentiateNormals >::SecondFundamentalFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals ) 
	: SecondFundamentalFormField( &vertices[0] , &normals[0] )
{}

template< unsigned int K , bool DifferentiateNormals >
typename SecondFundamentalFormField< K , DifferentiateNormals >::T
SecondFundamentalFormField< K , DifferentiateNormals >::Value
(
	const Point< double , Dim > & normal ,
	const Point< double , Dim > normals[K+1] ,
	const Matrix< double , K , Dim > & xForm ,
	Position< K > p
)
{
	if constexpr( DifferentiateNormals )
	{
		Differential< K , Point< double , Dim > > _dN = DNormalizedValue( LinearInterpolant< K , Point< double , Dim > >::Value( normals , p ) , LinearInterpolant< K , Point< double , Dim > >::DValue( normals , p ) );
		Matrix< double , K , Dim > dN;
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int n=0 ; n<K+1 ; n++ ) dN(k,n) = _dN[k][n];
		return PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::Value( normal , normals , xForm , p ).transpose() * dN;
	}
	else
	{
		T cov;
		Point< double , Dim > n = NormalizedValue( LinearInterpolant< K , Point< double , Dim > >::Value( normals , p ) );
		SimplexProcessing::Differential< K , Matrix< double , K , K+1 > > di2e = PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::DValue( normal , normals , xForm , p );

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
CovariantDerivativeField< K , N , VectorField >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const VectorField & VF ,
	Position< K > p
)
{
	SquareMatrix< double , K > t;

	// The transformation from extrinsic to intrinsic tangents
	Matrix< double , N , K > e2i = PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value( gInv , normal , normals , xForm , p );

	// The differentials of the extrinsic tangent vector fields
	Differential< K , Point< double , N > > dvf = VF.d(p);

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
CovariantDirectionalDerivativeField< K , N , DirectionField , VectorField >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const DirectionField & DF ,
	const VectorField & VF ,
	Position< K > p
)
{
	// The transformation from extrinsic to intrinsic tangents
	Matrix< double , N , K > e2i = PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value( gInv , normal , normals , xForm , p );

	// The evaluations of the vector fields, in the triangle tangent basis
	Point< double , K > dir = e2i( DF(p) );

	// The differentials of the extrinsic tangent vector fields
	SimplexProcessing::Differential< K , Point< double , N > > dvf = VF.d(p);

	// The covariant derivative
	return dvf( dir );
}

/////////////////////
// DivergenceField //
/////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
typename DivergenceField< K , N , VectorField >::T
DivergenceField< K , N , VectorField >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const VectorField & VF ,
	Position< K > p
)
{
	return CovariantDerivativeField< K , N , VectorField >::Value( gInv , normal , normals , xForm , VF , p ).trace();
}

////////////////////////////////////////
// CovariantDerivativeDifferenceField //
////////////////////////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField1 , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField2 >
typename CovariantDerivativeDifferenceField< K , N , VectorField1 , VectorField2 >::T
CovariantDerivativeDifferenceField< K , N , VectorField1 , VectorField2 >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const VectorField1 & VF1 ,
	const VectorField2 & VF2 ,
	Position< K > p
)
{
	// The transformation from extrinsic to intrinsic tangents
	Matrix< double , N , K > e2i = PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value( gInv , normal , normals , xForm , p );

	// The evaluations of the vector fields, in the triangle tangent basis
	Point< double , K > vf1 = e2i( VF1(p) ) , vf2 = e2i( VF2(p) );

	// The differentials of the extrinsic tangent vector fields
	Differential< K , Point< double , K+1 > > dvf1 = VF1.d(p) , dvf2 = VF2.d(p);

	// The difference of covariant derivatives
	return dvf2( vf1 ) - dvf1( vf2 );
}

/////////////////////
// LieBracketField //
/////////////////////
template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField1 , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField2 >
typename LieBracketField< K , N , VectorField1 , VectorField2 >::T
LieBracketField< K , N , VectorField1 , VectorField2 >::Value
(
	const SquareMatrix< double , K > & gInv ,
	const Point< double , N > & normal ,
	const Point< double , N > normals[K+1] ,
	const Matrix< double , K , N > & xForm ,
	const VectorField1 & VF1 ,
	const VectorField2 & VF2 ,
	Position< K > p
)
{
	// The intrinsic/extrinsic transformations
	Matrix< double , K , N > i2e = PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::Value( normal , normals , xForm , p );
	Matrix< double , N , K > e2i = PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::Value( gInv , i2e );
	Differential< K , Matrix< double , K , N > > di2e = PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::DValue( normal , normals , xForm , p );
	Differential< K , Matrix< double , N , K > > de2i = PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::DValue( gInv , di2e );

	// The extrinsic evaluations of the vector-fields and their derivatives
	Point< double , N > v1_e = VF1(p);
	Point< double , N > v2_e = VF2(p);
	Differential< K , Point< double , N > > dv1_e = VF1.d(p);
	Differential< K , Point< double , N > > dv2_e = VF2.d(p);

	// The intrinsic evaluations of the vector-fields and their derivatives
	Point< double , K > v1_i = IntrinsicVectorField< K , N , VectorField1 >::Value( e2i , v1_e );
	Point< double , K > v2_i = IntrinsicVectorField< K , N , VectorField2 >::Value( e2i , v2_e );
	Differential< K , Point< double , K > > dv1_i = IntrinsicVectorField< K , N , VectorField1 >::DValue( e2i , de2i , v1_e , dv1_e );
	Differential< K , Point< double , K > > dv2_i = IntrinsicVectorField< K , N , VectorField2 >::DValue( e2i , de2i , v2_e , dv2_e );

	// The intrinsic bracket
	Point< double , K > b_i;
	for( unsigned int k=0 ; k<K ; k++ ) b_i += dv2_i[k] * v1_i[k] - dv1_i[k] * v2_i[k];

	// The extrinsic bracket
	return i2e * b_i;
}