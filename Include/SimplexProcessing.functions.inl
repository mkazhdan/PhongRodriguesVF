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

template< unsigned int K >
template< typename T >
double InverseMetricTensor< K >::operator()( Differential< K , T > d1 , Differential< K , T > d2 ) const
{
	double d = 0;
	for( unsigned int k1=0 ; k1<K ; k1++ ) for( unsigned int k2=0 ; k2<K ; k2++ ) d += _gInv(k1,k2) * DotProduct( d1[k2] , d2[k1] );
	return d;
}

template< unsigned int K >
template< typename T >
Differential< K , T > InverseMetricTensor< K >::operator()( Differential< K , T > d ) const
{
	Differential< K , T > _d;
	for( unsigned int k1=0 ; k1<K ; k1++ )
	{
		T t{};
		for( unsigned int k2=0 ; k2<K ; k2++ ) t += _gInv(k2,k1) * d[k2];
		_d[k1] = t;
	}
	return _d;
}

template< unsigned int K , unsigned int D >
MetricTensor< K > MetricTensorFromEmbedding( const Point< double , D > v[K+1] )
{
	MetricTensor< K > g;
	for( unsigned int k1=0 ; k1<K ; k1++ ) for( unsigned int k2=0 ; k2<K ; k2++ ) g(k1,k2) = Point< double , D >::Dot( v[k1+1]-v[0] , v[k2+1]-v[0] );
	return g;
}

template< unsigned int D >
SquareMatrix< double , D > RodriguesRotation( Point< double , D > v , Point< double , D > w )
{
	// Normalize
	v /= Point< double , D >::Length( v );
	w /= Point< double , D >::Length( w );

	SquareMatrix< double , D > W = OuterProduct( w , v ) - OuterProduct( v , w );
	return SquareMatrix< double , D >::Identity() + W + W*W / ( 1. + Point< double , D >::Dot( v , w ) );
}

template< unsigned int D >
Differential< D , SquareMatrix< double , D > > D1RodriguesRotation( Point< double , D > v , Point< double , D > w )
{
	// R(v,w) = Id. + W(v,w) + ( W(v,w) * W(v,w) ) / ( 1 + <v,w> )
	//   with W(v,w) = w x v - v x w
	// 
	// Fixing w:
	//		W(v) = w x v - v x w
	//		dW(dN) = w x dN - dN x v
	// 
	//		W2(v) = ( w x v - v x w ) * ( w x v - v x w )
	//		dW2(dN) = [ w x dN - dN x w ] * W + W * [ w x dN - dN x w ]

	// Normalize
	double l = Point< double , D >::Length( v );
	v /= l;
	w /= Point< double , D >::Length( w );

	SquareMatrix< double , D > P = ( SquareMatrix< double , D >::Identity() - OuterProduct( v , v ) ) / l;

	SquareMatrix< double , D > W = OuterProduct( w , v ) - OuterProduct( v , w );
	Differential< D , SquareMatrix< double , D > > diff;
	for( unsigned int d=0 ; d<D ; d++ )
	{
		Point< double , D > dv;
		{
			Point< double , D > e;
			e[d] = 1;
			dv = P * e;
		}
		SquareMatrix< double , D > dW = OuterProduct( w , dv ) - OuterProduct( dv , w );
		double normalization = 1. + Point< double , D >::Dot( v , w );

		diff[d] = 
			// dW
			dW
			// d(W^2)/Normalization
			+ ( dW * W + W * dW ) / normalization
			// W^2*d(1/Normalization)
			- ( W * W * Point< double , D >::Dot( w , dv ) ) / ( normalization * normalization )
			;
	}

	return diff;
}

template< unsigned int D >
Differential< D , SquareMatrix< double , D > > D2RodriguesRotation( Point< double , D > v , Point< double , D > w )
{
	// R(v,w) = Id. + W(v,w) + ( W(v,w) * W(v,w) ) / ( 1 + <v,w> )
	//   with W(v,w) = w x v - v x w
	// 
	// Fixing v:
	//		W(w) = w x v - v x w
	//		dW(dN) = dN x v - v x dN
	// 
	//		W2(w) = ( w x v - v x w ) * ( w x v - v x w )
	//		dW2(dN) = [ dN x v - v x dN ] * W + W * [ dN x v - v x dN ]

	// Normalize
	double l = Point< double , D >::Length( w );
	v /= Point< double , D >::Length( v );
	w /= l;

	// The scaled projection accounting for normalization
	SquareMatrix< double , D > P = ( SquareMatrix< double , D >::Identity() - OuterProduct( w , w ) ) / l;

	SquareMatrix< double , D > W = OuterProduct( w , v ) - OuterProduct( v , w );
	double normalization = 1. + Point< double , D >::Dot( v , w );

	Differential< D , SquareMatrix< double , D > > diff;
	for( unsigned int d=0 ; d<D ; d++ )
	{
		Point< double , D > dw;
		{
			Point< double , D > e;
			e[d] = 1;
			dw = P * e;
		}
		SquareMatrix< double , D > dW = OuterProduct( dw , v ) - OuterProduct( v , dw );

		diff[d] = 
			// dW
			dW
			// d(W^2)/Normalization
			 + ( dW * W + W * dW ) / normalization
			// W^2*d(1/Normalization)
			- ( W * W * Point< double , D >::Dot( dw , v ) ) / ( normalization * normalization )
			;
	}

	return diff;
}

template< HasInnerProductSpace T >
double DotProduct( const T &v1 , const T &v2 )
{
	return T::Dot( v1 , v2 );
}

double DotProduct( double v1 , double v2 )
{
	return v1 * v2;
}

template< typename T , HasArray< T > Samples >
T SampleSum( Samples && S , unsigned int N )
{
	T v = {};
	for( unsigned int n=0 ; n<N ; n++ ) v += S[n];
	return v;
}

template< HasDotProduct T >
T NormalizedValue( const T & t )
{
	return t / sqrt( DotProduct( t , t ) );
}

template< unsigned int K , HasDotProduct T >
Differential< K , T > DNormalizedValue( T t , Differential< K , T > dt )
{
	// N(x) = F(x) / < F(x) , F(x) >^0.5
	// dN(x) = dF(x) / < F(x) , F(x) >^0.5 - 0.5 * F(x) / < F(x) , F(x) >^1.5 * 2 * dF(x)
	//       = dF(x) / < F(x) , F(x) >^0.5 - F(x) / < F(x) , F(x) >^1.5 * dF(x)

	double l = sqrt( DotProduct( t , t ) );
	t /= l;
	for( unsigned int k=0 ; k<K ; k++ )
	{
		// Normalize
		dt[k] /= l;
		// Project out the component in direction f
		dt[k] -= DotProduct( dt[k] , t ) * t;
	}
	return dt;
}

template< unsigned int K >
SquareMatrix< double , K > CovariantDerivative( AutoDiff::Tensor< K , K , K > gamma , Point< double , K > coords , Differential< K , Point< double , K > > dCoords )
{
	// Given a tangent vector field:
	//		v = \sum_k f_k * e_k,
	// the covariant derivative along direction e_i is
	//		dv/de_i = d(\sum_k f_k * e_k)/de_i
	//              = \sum_k df_k/de_i * e_k + f_k * de_k/de_i
	//              = \sum_k df_k/de_i * e_k + f_k * \sum_l e_l * \Gamma_{ki}^l

	SquareMatrix< double , K > d;

	for( unsigned int i=0 ; i<K ; i++ )
	{
		for( unsigned int k=0 ; k<K ; k++ ) d(i,k) += dCoords[i][k];
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int l=0 ; l<K ; l++ ) d(i,l) += coords[k] * gamma(k,i,l);
	}
	return d;
}

template< unsigned int K , HasSimplexFunction< K , ScaleFactor > ScaleFactorField , HasSimplexInvocable< K > Field >
auto ScaledField( ScaleFactorField && SF , Field && F )
{
	return [_SF=std::forward< ScaleFactorField >(SF),_F=std::forward< Field >(F)]( Position< K > p ){ return _SF(p) * _F(p); };
}

template< unsigned int K >
double UnitWeightField< K >::operator()( Position< K > ) const { return 1.; }

template< unsigned int K , unsigned int N >
SquareMatrix< double , N > IdentityField< K , N >::operator()( Position< K > ) const { return SquareMatrix< double , N >::Identity(); };
