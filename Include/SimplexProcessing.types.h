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

#ifndef SIMPLEX_PROCESSING_TYPES_INCLUDED
#define SIMPLEX_PROCESSING_TYPES_INCLUDED

namespace MishaK
{
	namespace SimplexProcessing
	{
		using Measure = double;
		using ScaleFactor = double; // Used to transform one measure to another (via multiplication)
		using SimplexIndex = unsigned int;

		template< unsigned int K >                    using Position = Point< double , K >;
		template< unsigned int K , typename T >       using Differential = Point< T , K , double >;
		template< unsigned int K >                    using MetricTensor = SquareMatrix< double , K >;
		template< unsigned int K , unsigned int Dim > using EmbeddedSimplex = Simplex< double , Dim , K >;
		template< unsigned int K >
		struct InverseMetricTensor /* : public BilinearForm< T > */
		{
			InverseMetricTensor( void ) : _gInv( SquareMatrix< double , K >::Identity() ){}
			InverseMetricTensor( MetricTensor< K > g ) : _gInv( g.inverse() ){}
			template< typename T >
			double operator()( Differential< K , T > d1 , Differential< K , T > d2 ) const;
			template< typename T >
			Differential< K , T > operator()( Differential< K , T > d ) const;
		protected:
			SquareMatrix< double , K > _gInv;
		};

		// Takes a class with a function operator taking an integer index and wraps it to look like an array.
		// [NOTE] The function captures by value.
		template< typename Accessor >
		struct ArrayWrapper
		{
			ArrayWrapper( const Accessor & accessor ) : _accessor( accessor ){}
			decltype( std::declval< Accessor >()( std::declval< size_t >() ) ) operator[]( size_t idx ) const { return _accessor(idx); }
		protected:
			Accessor _accessor;
		};
	}
}
#endif // SIMPLEX_PROCESSING_TYPES_INCLUDED
