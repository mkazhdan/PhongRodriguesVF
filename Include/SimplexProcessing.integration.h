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

#ifndef SIMPLEX_PROCESSING_INTEGRATOR_INCLUDED
#define SIMPLEX_PROCESSING_INTEGRATOR_INCLUDED

namespace MishaK
{
	namespace SimplexProcessing
	{
		namespace SystemIntegration
		{
			// Functionality for constructing system constraints, expressed relative to the test functions
			template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexSystemVectorFunction< K , NumF , T > SystemVectors >
			auto SystemVectorField( TestFunctions && Fs , SystemVectors && S );

			template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexlinearMapFunction< K , T > LinearMaps >
			auto SystemVectorField( Function && F , TestFunctions && Fs , LinearMaps && L );

			template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexBilinearFormFunction< K , T > BilinearForms >
				requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > )
			auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B );

			// Functionality for constructing system constraints, expressed relative to the test functions
			template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexSystemMatrixFunction< K , NumF , T > SystemMatrices >
				requires( !HasSimplexlinearMapFunction< SystemMatrices , K , T > && !HasSimplexBilinearFormFunction< SystemMatrices , K , T > )
			auto SystemMatrixField( TestFunctions && Fs , SystemMatrices && S );

			template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexlinearMapFunction< K , T > LinearMaps >
			auto SystemMatrixField( TestFunctions && Fs , LinearMaps && L );

			template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexBilinearFormFunction< K , T > BilinearForms >
				requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > )
			auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B );

			// A dimension-templated wrapper for the segment/triangle/tet integrators
			template< unsigned int K , unsigned int QuadratureSamples >
			struct MCIntegrator /* : public Sample< Measure > */
			{
				static SimplexProcessing::Position< K > Position( unsigned int n ){ return _Integrator::Positions[n]; }
				static SimplexProcessing::Measure Measure( unsigned int n ){ return _Integrator::Weights[n] / Factorial< K >(); }

				SimplexProcessing::Measure operator()( unsigned int n ) const { return Measure(n); }

				// The integral of a function over the simplex
				template< typename T , HasSimplexFunction< K , T > Function >
				static T Integral( Function && F );

			protected:
				static_assert( K==1 || K==2 || K==3 , "[ERROR] Only 1- , 2-, and 3-dimensional simplices supported" );
				using _Integrator = SimplexIntegrator< K , QuadratureSamples >;
			};
#include "SimplexProcessing.integration.inl"
		}
	}
}
#endif // SIMPLEX_PROCESSING_INTEGRATOR_INCLUDED
