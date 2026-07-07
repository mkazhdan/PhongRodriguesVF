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

#ifndef SIMPLEX_PROCESSING_FUNCTIONS_INCLUDED
#define SIMPLEX_PROCESSING_FUNCTIONS_INCLUDED

namespace MishaK
{
	namespace SimplexProcessing
	{
		// Returns the factorial of I
		template< unsigned int I >
		static constexpr unsigned int Factorial( void )
		{
			if constexpr( I==0 ) return 1;
			else return Factorial< I-1 >() * I;
		}

		// Returns the center of the simplex
		template< unsigned int K >
		Position< K > CenterPosition( void ){ Position< K > p ; for( unsigned int k=0 ; k<K ; k++ ) p[k] = 1./(K+1) ; return p; }

		// Computes the scale factor from the metric tensor
		template< unsigned int K >
		ScaleFactor ScaleFactorFromMetricTensor( MetricTensor< K > g ){ return sqrt( g.determinant() ); }

		// Computes the metric tensor from a simplex embedding
		template< unsigned int K , unsigned int D >
		MetricTensor< K > MetricTensorFromEmbedding( const Point< double , D > v[K+1] );

		// Computes the smallest rotation taking v to w
		template< unsigned int D >
		SquareMatrix< double , D > RodriguesRotation( Point< double , D > v , Point< double , D > w );

		// Computes the differential of the smallest rotation taking v to w, with respect to the first argument
		template< unsigned int D >
		Differential< D , SquareMatrix< double , D > > D1RodriguesRotation( Point< double , D > v , Point< double , D > w );

		// Computes the differential of the smallest rotation taking v to w, with respect to the second argument
		template< unsigned int D >
		Differential< D , SquareMatrix< double , D > > D2RodriguesRotation( Point< double , D > v , Point< double , D > w );

		// A general function return the dot-product between two elements of type T, with T being an inner-product space
		template< HasInnerProductSpace T >
		double DotProduct( const T & v1 , const T & v2 );

		// A specialization of the dot-product function for inputs that are doubles
		double DotProduct( double v1 , double v2 );

		// The identity map
		template< typename T > T Identity( const T & v ){ return v; }

		// Basic functionality for summing samples with a variable number of samples
		template< typename T , HasArray<T > Samples >
		T SampleSum( Samples && S , unsigned int N );

		template< HasDotProduct T >
		T NormalizedValue( const T & t );

		template< unsigned int K , HasDotProduct T >
		Differential< K , T > DNormalizedValue( T t , Differential< K , T > dt );

		// Computes the covariant derivative using the connection coefficients, the coefficients, and their derivatives
		template< unsigned int K >
		SquareMatrix< double , K > CovariantDerivative( AutoDiff::Tensor< K , K , K > gamma , Point< double , K > coords , Differential< K , Point< double , K > > dCoords );

		// Takes a scale-factor field and a field, and return the scaled field
		template< unsigned int K , HasSimplexFunction< K , ScaleFactor > ScaleFactorField , HasSimplexInvocable< K > Field >
		auto ScaledField( ScaleFactorField && SF , Field && F );

#include "SimplexProcessing.functions.inl"
	}
}
#endif // SIMPLEX_PROCESSING_FUNCTIONS_INCLUDED
