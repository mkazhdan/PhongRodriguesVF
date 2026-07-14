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

////////////////////////////
// ScalarSystem::Elements //
////////////////////////////
template< unsigned int K >
ScalarSystem< K >::Elements::Elements( void )
{
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		double x[K+1] = {0};
		x[k] = 1;
		_f[k] = LinearInterpolant< K , double >( x );
	}
}

//////////////////
// ScalarSystem //
//////////////////
template< unsigned int K >
template< unsigned int N , HasSimplexFunction< K , Differential< K , double > > VectorField , bool MeasureScale >
auto ScalarSystem< K >::DerivationMatrix( const Point< double , N > vertices[K+1] , const VectorField & VF )
{
	return [ elements=Elements() , measureScale=MeasureScaleField<K>( vertices ) , vf=VF ]( Position< K > p )
		{
			double values[NumElements];
			Differential< K , double > dValues[NumElements];

			Differential< K , double > v = vf(p);

			for( unsigned int e=0 ; e<NumElements ; e++ )
			{
				auto element = elements[e];
				values[e] = element(p);
				dValues[e] = element.d(p);
			}
			SquareMatrix< double , NumElements > D;
			for( unsigned int m=0 ; m<NumElements ; m++ )
			{
				// Compute the derivative of the m-th element function along the vector-field
				double value = 0;
				for( unsigned int k=0 ; k<K ; k++ ) value += dValues[m][k] * v[k];

				// Integrate that against the n-the element function
				for( unsigned int n=0 ; n<NumElements ; n++ ) D(n,m) = value * values[n];
			}

			if constexpr( MeasureScale ) return D * measureScale( p );
			else                         return D;
		};
}

template< unsigned int K >
template< unsigned int N , HasSimplexFunction< K , double > Function , bool MeasureScale >
auto ScalarSystem< K >::MassVector( const Point< double , N > vertices[K+1] , const Function & F )
{
	return [ elements=Elements() , measureScale=MeasureScaleField<K>( vertices ) , f=F ]( Position< K > p )
		{
			Point< double , NumElements > m;
			double value = f(p);
			for( unsigned int e=0 ; e<NumElements ; e++ ) m[e] = value * elements[e](p);
			if constexpr( MeasureScale ) return m * measureScale( p );
			else                         return m;
		};
}  

template< unsigned int K >
template< unsigned int N , HasSimplexDifferentiableFunction< K , double > DifferentiableFunction , bool MeasureScale >
auto ScalarSystem< K >::StiffnessVector( const Point< double , N > vertices[K+1] , const DifferentiableFunction & F )
{
	return [ elements=Elements() , gInv=std::conditional_t< MeasureScale , MeasureScaledInverseFirstFundamentalFormField< K > , InverseFirstFundamentalFormField< K > >( vertices ) , f=F ]( Position< K > p )
		{
			Point< double , NumElements > s;
			Differential< K , double > dValue = gInv( p ) * f.d(p);
			for( unsigned int e=0 ; e<NumElements ; e++ ) s[e] = DotProduct( dValue , elements[e].d(p) );
			return s;
		};
}

template< unsigned int K >
template< unsigned int N , HasSimplexFunction< K , Differential< K , double > > VectorField , bool MeasureScale >
auto ScalarSystem< K >::StiffnessVector( const Point< double , N > vertices[K+1] , const VectorField & VF )
{
	return [ elements=Elements() , gInv=std::conditional_t< MeasureScale , MeasureScaledInverseFirstFundamentalFormField< K > , InverseFirstFundamentalFormField< K > >( vertices ) , vf=VF ]( Position< K > p )
		{
			Point< double , NumElements > s;
			Differential< K , double > dValue = gInv( p ) * vf(p);
			for( unsigned int e=0 ; e<NumElements ; e++ ) s[e] = DotProduct( dValue , elements[e].d(p) );
			return s;
		};
}

template< unsigned int K >
template< unsigned int N , bool MeasureScale >
auto ScalarSystem< K >::MassMatrix( const Point< double , N > vertices[K+1] )
{
	return [ elements=Elements() , measureScale=MeasureScaleField<K>( vertices ) ]( Position< K > p )
		{
			double values[NumElements];
			for( unsigned int e=0 ; e<NumElements ; e++ ) values[e] = elements[e]( p );
			SquareMatrix< double , NumElements > M;
			for( unsigned int i=0 ; i<NumElements ; i++ ) for( unsigned int j=0 ; j<NumElements ; j++ ) M(i,j) = values[i] * values[j];
			if constexpr( MeasureScale ) return M * measureScale( p );
			else                         return M;
		};
}

template< unsigned int K >
template< unsigned int N , bool MeasureScale >
auto ScalarSystem< K >::StiffnessMatrix( const Point< double , N > vertices[K+1] )
{
	return [ elements=Elements() , GInv=std::conditional_t< MeasureScale , MeasureScaledInverseFirstFundamentalFormField< K > , InverseFirstFundamentalFormField< K > >(vertices) ]( Position< K > p )
		{
			Differential< K , double > dValues[NumElements];
			SquareMatrix< double , K > gInv = GInv( p );
			for( unsigned int e=0 ; e<NumElements ; e++ ) dValues[e] = elements[e].d( p );
			SquareMatrix< double , NumElements > S;
			for( unsigned int i=0 ; i<NumElements ; i++ )
			{
				Differential< K , double > dValue = gInv( dValues[i] );
				for( unsigned int j=0 ; j<NumElements ; j++ ) S(i,j) = DotProduct( dValue , dValues[j] );
			}
			return S;
		};
}


////////////////////////////////////
// PhongRodriguesSystem::Elements //
////////////////////////////////////
template< unsigned int K , unsigned int N >
PhongRodriguesSystem< K , N >::Elements::Elements( const Point< double , N > normals[K+1] )
{
	Point< double , N > vf[K+1];
	for( unsigned int e=0 , k=0 ; k<=K ; k++ ) for( unsigned int n=0 ; n<N ; n++ , e++ )
	{
		vf[k][n] = 1;
		_f[e] = PhongRodriguesVectorField< K , N >( normals , vf );
		vf[k][n] = 0;
	}
}

//////////////////////////
// PhongRodriguesSystem //
//////////////////////////
template< unsigned int K , unsigned int N >
template< HasSimplexFunction< K , Point< double , N > > VectorField , bool MeasureScale , HasSimplexFunction< K , SquareMatrix< double , N > > InnerProductField >
auto PhongRodriguesSystem< K , N >::MassVector( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & VF , InnerProductField && IP )
{
	return [ elements=Elements( normals ) , measureScale=MeasureScaleField<K>( vertices ) , vf=VF , ip=IP ]( Position< K > p )
		{
			Point< double , NumElements > m;
			Point< double , N > v = vf(p);
			if constexpr( !std::same_as< InnerProductField , IdentityField< K , N > > ) v = ip(p) * v;
			for( unsigned int e=0 ; e<NumElements ; e++ ) m[e] = Point< double , N >::Dot( v , elements[e](p) );
			if constexpr( MeasureScale ) return m * measureScale( p );
			else                         return m;
		};
}

template< unsigned int K , unsigned int N >
template< bool MeasureScale , HasSimplexFunction< K , SquareMatrix< double , N > > InnerProductField >
auto PhongRodriguesSystem< K , N >::MassMatrix( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , InnerProductField && IP )
{
	return [ elements=Elements( normals ) , measureScale=MeasureScaleField<K>( vertices ) , ip=IP ]( Position< K > p )
		{
			Point< double , N > values[NumElements];
			for( unsigned int e=0 ; e<NumElements ; e++ ) values[e] = elements[e]( p );
			SquareMatrix< double , NumElements > M;
			if constexpr( std::same_as< InnerProductField , IdentityField< K , N > > )
			{
				for( unsigned int i=0 ; i<NumElements ; i++ ) for( unsigned int j=0 ; j<NumElements ; j++ )
					M(i,j) = DotProduct( values[i] , values[j] );
			}
			else
			{
				SquareMatrix< double , N > I = ip(p);
				for( unsigned int i=0 ; i<NumElements ; i++ )
				{
					Point< double , N > value = I * values[i];
					for( unsigned int j=0 ; j<NumElements ; j++ ) M(i,j) = DotProduct( value , values[j] );
				}
			}
			if constexpr( MeasureScale ) return M * measureScale( p );
			else                         return M;
		};
}

template< unsigned int K , unsigned int N >
template< bool MeasureScale >
auto PhongRodriguesSystem< K , N >::StiffnessMatrix( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] )
{
	return [ elements=Elements( normals ) , GInv=std::conditional_t< MeasureScale , MeasureScaledInverseFirstFundamentalFormField< K > , InverseFirstFundamentalFormField< K > >( vertices ) , normalField=NormalizationField< K , Point< double , N > , LinearInterpolant< K , Point< double , N > > >( normals ) ]( Position< K > p )
		{
			Matrix< double , K , N > covariantDerivatives[ PhongRodriguesSystem< K ,N >::NumElements ];
			Point< double , N > normal = normalField(p);
			SquareMatrix< double , K > gInv = GInv( p );

			SquareMatrix< double , PhongRodriguesSystem< K ,N >::NumElements > S;

			for( unsigned int e=0 ; e<PhongRodriguesSystem< K ,N >::NumElements ; e++ )
			{
				Differential< K , Point< double , N > > dValue = elements[e].d(p);
				for( unsigned int k=0 ; k<K ; k++ )
				{
					dValue[k] -= normal * Point< double , N >::Dot( dValue[k] , normal );
					for( unsigned int n=0 ; n<N ; n++ ) covariantDerivatives[e](k,n) = dValue[k][n];
				}
			}
			for( unsigned int n=0 ; n<PhongRodriguesSystem< K ,N >::NumElements ; n++ ) for( unsigned int m=0 ; m<PhongRodriguesSystem< K ,N >::NumElements ; m++ )
				S(n,m) = SquareMatrix< double , K >::Trace( gInv * covariantDerivatives[n].transpose() * covariantDerivatives[m] );
			return S;
		};
}

template< unsigned int K , unsigned int N >
template< unsigned int Components , bool MeasureScale >
auto PhongRodriguesSystem< K , N >::ComponentStiffnessMatrix( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] )
{
	return [ elements=Elements( normals ) , I2E = PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >( vertices , normals ) , G = FirstFundamentalFormField< K >( vertices ) , GInv=std::conditional_t< MeasureScale , MeasureScaledInverseFirstFundamentalFormField< K > , InverseFirstFundamentalFormField< K > >( vertices ) , measureScale=MeasureScaleField< K >( vertices ) ]( Position< K > p )
		{
			SquareMatrix< double , K > covariantDerivativeForms[ PhongRodriguesSystem< K ,N >::NumElements ];
			Matrix< double , N , K > i2e_t = I2E( p ).transpose();
			SquareMatrix< double , K > g = G( p );
			SquareMatrix< double , K > gInv = GInv( p );
			if constexpr( MeasureScale ) g /= measureScale(p);

			SquareMatrix< double , PhongRodriguesSystem< K ,N >::NumElements > S;

			for( unsigned int e=0 ; e<PhongRodriguesSystem< K ,N >::NumElements ; e++ )
			{
				SquareMatrix< double , K > form;
				Differential< K , Point< double , N > > dValue = elements[e].d(p);
				for( unsigned int k=0 ; k<K ; k++ )
				{
					Point< double , K > dual = i2e_t * dValue[k];
					for( unsigned int l=0 ; l<K ; l++ ) form(k,l)= dual[l];
				}

				if      constexpr( Components & StiffnessComponent::Curl()                                                             ) covariantDerivativeForms[e] += ( form - form.transpose() ) / 2.;
				if      constexpr( Components & StiffnessComponent::Divergence() && Components & StiffnessComponent::AntiHolomorphic() ) covariantDerivativeForms[e] += ( form + form.transpose() ) / 2.;
				else if constexpr( Components & StiffnessComponent::Divergence()                                                       ) covariantDerivativeForms[e] += g * ( gInv * form ).trace() / 2.;
				else if constexpr( Components & StiffnessComponent::AntiHolomorphic()                                                  ) covariantDerivativeForms[e] += ( form + form.transpose() ) / 2. - g * ( gInv * form ).trace() / 2.;
			}
			// The inner product of two maps M,N: TM -> T^*M is
			//		<M,N> = tr( ( gInv * ( gInv * N ).transpose() * g * ( gInv * M ) );
			//            = tr( ( gInv * N.transpose() * gInv * M );
			for( unsigned int n=0 ; n<PhongRodriguesSystem< K ,N >::NumElements ; n++ ) for( unsigned int m=0 ; m<PhongRodriguesSystem< K ,N >::NumElements ; m++ )
				S(n,m) = SquareMatrix< double , K >::Trace( gInv * covariantDerivativeForms[n].transpose() * gInv * covariantDerivativeForms[m] );

			if constexpr( MeasureScale ) return S / measureScale( p );
			else                         return S;
		};
}


