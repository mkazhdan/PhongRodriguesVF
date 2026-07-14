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

#ifndef SIMPLEX_PROCESSING_CONCEPTS_INCLUDED
#define SIMPLEX_PROCESSING_CONCEPTS_INCLUDED

namespace MishaK
{
	namespace SimplexProcessing
	{
		template< typename T >
		concept HasInnerProductSpace = std::derived_from< T , InnerProductSpace< double , T > >;

		template< typename T >
		concept HasDotProduct = std::same_as< T , double > || HasInnerProductSpace< T >;

		template< typename Array , typename T >
		concept HasArray = requires( const Array array , size_t idx ) { { array[idx] } -> std::same_as< T >; };

		template< typename Function , typename ... Args >
		concept HasInvocable = std::invocable< Function , Args... >;

		template< typename Function , unsigned int K >
		concept HasSimplexInvocable = HasInvocable< Function , Position< K > >;

		template< typename Function , typename RetValue , typename ... Args >
		concept HasFunction = requires( const Function f , Args ... args ) { { f(args...) } -> std::same_as< RetValue >; };

		template< typename Field , unsigned int K , typename T >
		concept HasSimplexFunction = HasFunction< Field , T , Position< K > >;

		template< typename Field , unsigned int K , typename T >
		concept HasSimplexDifferentiableFunction = HasSimplexFunction< Field , K , T > && requires( const Field f , Position< K > p ) { { f.d(p) } -> std::same_as< Differential< K , T > >; };

		template< typename Field , unsigned int K >
		concept HasArrayOfSimplexInvocables = requires( const Field f , size_t idx , Position< K > p ) { { f[idx](p) }; };

		template< typename Field , unsigned int K , typename T >
		concept HasArrayOfSimplexFunctions = requires( const Field f , size_t idx , Position< K > p ) { { f[idx](p) } -> std::same_as< T >; };

		template< typename Field , unsigned int K , typename T >
		concept HasArrayOfSimplexDifferentiableFunctions = requires( const Field f , size_t idx , Position< K > p ) { { f[idx].d(p) } -> std::same_as< Differential< K , T > >; };
	}
}
#endif // SIMPLEX_PROCESSING_CONCEPTS_INCLUDED
