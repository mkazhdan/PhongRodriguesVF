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
		/////////////////////////////////////
		// Covariant derivative components //
		/////////////////////////////////////
		struct StiffnessComponent
		{
			static constexpr unsigned int Divergence     ( void ){ return _Divergence                           ; }
			static constexpr unsigned int Curl           ( void ){ return               _Curl                   ; }
			static constexpr unsigned int AntiHolomorphic( void ){ return                       _AntiHolomorphic; }
			static constexpr unsigned int Connection     ( void ){ return _Divergence | _Curl | _AntiHolomorphic; }
			static constexpr unsigned int Holomorphic    ( void ){ return _Divergence | _Curl                   ; }
			static constexpr unsigned int Hodge          ( void ){ return _Divergence | _Curl                   ; }
			static constexpr unsigned int Killing        ( void ){ return _Divergence |         _AntiHolomorphic; }

			inline static const std::vector< std::string > Names = { "divergence" , "curl" , "anti-holomorphic" , "connection" , "holomorphic" , "hodge" , "killing" };
		protected:
			enum
			{
				_Divergence = 1 ,
				_Curl = 2 ,
				_AntiHolomorphic = 4
			};
		};

		///////////////////////////////////////////////
		// Scalar system matrix/vector functionality //
		///////////////////////////////////////////////
		template< unsigned int K >
		struct ScalarSystem
		{
			static const unsigned int NumElements = K+1;

			///////////////////////////////////
			// First-order Lagrange elements //
			///////////////////////////////////
			struct Elements
			{
				Elements( void );
				LinearInterpolant< K , double > operator[]( size_t k ) const { return _f[k]; }
			protected:
				LinearInterpolant< K , double > _f[NumElements];
			};

			//////////////////////////
			// Vector functionality //
			//////////////////////////

			// Returns the functionality for computing the dual vector defined by integrating the function against the linear elements
			template< unsigned int N , HasSimplexFunction< K , double > Function , bool MeasureScale=true >
			static auto MassVector( const Point< double , N > vertices[K+1] , const Function & F );

			template< unsigned int N , HasSimplexFunction< K , double > Function , bool MeasureScale=true >
			static auto MassVector( const Simplex< double , N , K > & vertices , const Function & F ){ return MassVector< N , Function , MeasureScale >( &vertices[0] , F ); }

			// Returns the functionality for computing the dual vector defined by integrating the differential function against the differentials of the linear elements
			template< unsigned int N , HasSimplexDifferentiableFunction<  K , double > DifferentiableFunction , bool MeasureScale=true >
			static auto StiffnessVector( const Point< double , N > vertices[K+1] , const DifferentiableFunction & F );

			template< unsigned int N , HasSimplexDifferentiableFunction< K , double > DifferentiableFunction , bool MeasureScale=true >
			static auto StiffnessVector( const Simplex< double , N , K > & vertices , const DifferentiableFunction & F ){ return StiffnessVector< N , DifferentiableFunction , MeasureScale >( &vertices[0] , F ); }

			// Returns the functionality for computing the dual vector defined by integrating the vector-field function against the differentials of the linear elements
			template< unsigned int N , HasSimplexFunction< K , Differential< K , double > > VectorField , bool MeasureScale=true >
			static auto StiffnessVector( const Point< double , N > vertices[K+1] , const VectorField & VF );

			template< unsigned int N , HasSimplexFunction< K , Differential< K , double > > VectorField , bool MeasureScale=true >
			static auto StiffnessVector( const Simplex< double , N , K > & vertices , const VectorField & VF ){ return StiffnessVector< N , VectorField , MeasureScale >( &vertices[0] , VF ); }

			//////////////////////////
			// Matrix functionality //
			//////////////////////////

			// Returns the functionality for computing the derivation matrix associated to an intrinsic tangent vector-field
			template< unsigned int N , HasSimplexFunction< K , Differential< K , double > > VectorField , bool MeasureScale=true >
			static auto DerivationMatrix( const Point< double , N > vertices[K+1] , const VectorField & VF );

			template< unsigned int N , HasSimplexFunction< K , Differential< K , double > > VectorField , bool MeasureScale=true >
			static auto DerivationMatrix( const Simplex< double , N , K > & vertices , const VectorField & VF ){ return DerivationMatrix< N , VectorField , MeasureScale >( &vertices[0] , VF ); }

			// Returns the functionality for computing the mass matrix
			template< unsigned int N , bool MeasureScale=true >
			static auto MassMatrix( const Point< double , N > vertices[K+1] );

			template< unsigned int N , bool MeasureScale=true >
			static auto MassMatrix( const Simplex< double , N , K > & vertices ){ return MassMatrix< N , MeasureScale >( &vertices[0] ); }

			// Returns the functionality for computing the stiffness matrix
			template< unsigned int N , bool MeasureScale=true >
			static auto StiffnessMatrix( const Point< double , N > vertices[K+1] );

			template< unsigned int N , bool MeasureScale=true >
			static auto StiffnessMatrix( const Simplex< double , N , K > & vertices ){ return StiffnessMatrix< N , MeasureScale >( &vertices[0] ); }
		};

		template< unsigned int K , unsigned int N >
		struct PhongRodriguesSystem
		{
			static const unsigned int NumElements = (K+1)*N;

			/////////////////////////////////////////////////////
			// Phong-Rodrigues extrinsic vector-field elements //
			/////////////////////////////////////////////////////
			template< bool Modulate >
			struct _Elements
			{
				_Elements( const Point< double , N > n[K+1] );
				PhongRodriguesVectorField< K , N , Modulate > operator[]( size_t k ) const { return _f[k]; }
			protected:
				PhongRodriguesVectorField< K , N , Modulate > _f[NumElements];
			};

			using Elements = _Elements< true >;
			using UnomdulatedElements = _Elements< false >;

			//////////////////////////
			// Vector functionality //
			//////////////////////////

			// Returns the functionality for computing the dual vector defined by integrating the function against the Phong-Rodrigues elements
			template< HasSimplexFunction< K , Point< double , N > > VectorField , bool MeasureScale=true , HasSimplexFunction< K , SquareMatrix< double , N > > InnerProductField=IdentityField< K , N > >
			static auto MassVector( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & VF , InnerProductField && IP=IdentityField< K , N >() );

			template< HasSimplexFunction< K , Point< double , N > > VectorField , bool MeasureScale=true , HasSimplexFunction< K , SquareMatrix< double , N > > InnerProductField=IdentityField< K , N > >
			static auto MassVector( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & VF , InnerProductField && IP=IdentityField< K , N >() ){ return MassVector< VectorField , MeasureScale >( &vertices[0] , &normals[0] , VF , std::forward< InnerProductField >( IP ) ); }

			// Returns the functionality for computing the mass matrix
			template< bool MeasureScale=true , HasSimplexFunction< K , SquareMatrix< double , N > > InnerProductField=IdentityField< K , N > >
			static auto MassMatrix( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , InnerProductField && IP=IdentityField< K , N >() );

			template< bool MeasureScale=true , HasSimplexFunction< K , SquareMatrix< double , N > > InnerProductField=IdentityField< K , N > >
			static auto MassMatrix( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , InnerProductField && IP=IdentityField< K , N >() ){ return MassMatrix< MeasureScale , InnerProductField >( &vertices[0] , &normals[0] , std::forward< InnerProductField >( IP ) ); }

			// Returns the functionality for computing the (covariant) stiffness matrix
			template< bool MeasureScale=true >
			static auto StiffnessMatrix( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] );

			template< bool MeasureScale=true >
			static auto StiffnessMatrix( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals ){ return StiffnessMatrix< MeasureScale >( &vertices[0] , &normals[0] ); }

			// Returns the functionality for computing the components of the covariant stiffness matrix
			template< unsigned int Components , bool MeasureScale=true >
			static auto ComponentStiffnessMatrix( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] );

			template< unsigned int Components , bool MeasureScale=true >
			static auto ComponentStiffnessMatrix( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals ){ return ComponentStiffnessMatrix< Components , MeasureScale >( &vertices[0] , &normals[0] ); }
		};

#include "SimplexProcessing.elements.inl"
	}
}
#endif // SIMPLEX_PROCESSING_ELEMENTS_INCLUDED
