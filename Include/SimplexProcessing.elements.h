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

#ifndef SIMPLEX_PROCESSING_ELEMENTS_INCLUDED
#define SIMPLEX_PROCESSING_ELEMENTS_INCLUDED

namespace MishaK
{
	namespace SimplexProcessing
	{
		template< unsigned int K >
		struct LinearElements
		{
			static const unsigned int NumElements = K+1;

			LinearElements( void );
			LinearInterpolant< K , double > operator[]( size_t k ) const { return _f[k]; }
		protected:
			LinearInterpolant< K , double > _f[NumElements];
		};

		template< unsigned int K , unsigned int N , bool Modulate=true >
		struct PhongRodriguesVectorElements
		{
			static const unsigned int NumElements = (K+1)*N;

			PhongRodriguesVectorElements( const Point< double , N > n[K+1] );
			PhongRodriguesVectorField< K , N , Modulate > operator[]( size_t k ) const { return _f[k]; }
		protected:
			PhongRodriguesVectorField< K , N , Modulate > _f[NumElements];
		};

#include "SimplexProcessing.elements.inl"
	}
}
#endif // SIMPLEX_PROCESSING_ELEMENTS_INCLUDED
