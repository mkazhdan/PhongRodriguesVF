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

#ifndef SIMPLICIAL_MESH_PROCESSING_FIELDS_INCLUDED
#define SIMPLICIAL_MESH_PROCESSING_FIELDS_INCLUDED

namespace MishaK
{
	namespace SimplicialMesh
	{
		////////////////////////////////////////////////////////////////////////////////
		// Classes for wrapping per-simplex functionality into per-mesh functionality //
		////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , bool NeedsVertices , bool NeedsNormals , typename SimplexField , typename ... InTypes >
		struct _MeshField
		{
			struct Empty
			{
				Empty( void ){}
				template< typename T > Empty( const T & ){}
			};

			static const unsigned int Dim = K+1;
			using T = SimplexField::T;

			_MeshField( const EmbeddedMesh< K , Dim > & mesh , const std::vector< InTypes > & ... in ) requires( !NeedsNormals );
			_MeshField( const EmbeddedPhongMesh< K > & mesh , const std::vector< InTypes > & ... in );
			auto operator[]( size_t sIdx ) const;

		protected:
			const std::vector< SimplexIndex< K > > & _simplices;
			std::conditional_t< NeedsVertices , const std::vector< Point< double , Dim > > & , Empty > _vertices;
			std::conditional_t< NeedsNormals , const std::vector< Point< double , Dim > > & , Empty > _normals;
			std::tuple< const std::vector< InTypes > & ... > _in;
		};

		template< unsigned int K , bool NeedsVertices , bool NeedsNormals , typename SimplexField , typename ... InTypes >
		struct _DifferentiableMeshField : public _MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... >
		{
			// The differential function type
			struct Differential
			{
				Differential( const _MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... > & f ) : _f(f){}
				auto operator[]( size_t sIdx ) const { return _f[sIdx].d(); }
			protected:
				_MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... > _f;
			};

			// Inherit base constructors
			using _MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... >::_MeshField;

			// Return a function evaluating the differential
			Differential d( void ) const { return Differential( static_cast< const _MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... > & >( *this ) ); }
		};

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using MeshField = _MeshField< K , false , false , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using EmbeddedMeshField = _MeshField< K , true , false , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using PhongMeshField = _MeshField< K , false , true , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using EmbeddedPhongMeshField = _MeshField< K , true , true , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using DifferentiableMeshField = _DifferentiableMeshField< K , false , false , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using DifferentiableEmbeddedMeshField = _DifferentiableMeshField< K , true , false , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using DifferentiablePhongMeshField = _DifferentiableMeshField< K , false , true , SimplexField , InTypes ... >;

		template< unsigned int K , typename SimplexField , typename ... InTypes >
		using DifferentiableEmbeddedPhongMeshField = _DifferentiableMeshField< K , true , true , SimplexField , InTypes ... >;

		////////////////////////////////////////////////////////////////
		// A class returning the piecewise-linear interpolation field //
		////////////////////////////////////////////////////////////////
		template< unsigned int K , typename T=double >
		using ScalarField = DifferentiableMeshField< K , SimplexProcessing::LinearInterpolant< K , T > , T >;

		///////////////////////////////////////////////////
		// A class returning the second-fundamental form //
		///////////////////////////////////////////////////
		template< unsigned int K , bool DifferentiateNormals=true >
		using SecondFundamentalFormField = EmbeddedPhongMeshField< K , SimplexProcessing::SecondFundamentalFormField< K , DifferentiateNormals > >;

		///////////////////////////////////////////////
		// A class returning connection coefficients //
		///////////////////////////////////////////////
		template< unsigned int K >
		using ConnectionField = EmbeddedPhongMeshField< K , SimplexProcessing::ConnectionCoefficientField< K > >;

		////////////////////////////////////////////////////////////////////////////////////////////
		// Classes returning the transformation field between intrinsic/extrensic tangent vectors //
		////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesExtrinsicToIntrinsicTangentXFormField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > >;

		template< unsigned int K >
		using PhongRodriguesIntrinsicToExtrinsicTangentXFormField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K > >;

		///////////////////////////////////////////////////////////////////
		// A class returning the associated Phong-Rodrigues vector field //
		///////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesVectorField = PhongMeshField< K , SimplexProcessing::PhongRodriguesVectorField< K , K+1 > , Point< double , K+1 > >;


		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		// A class returning the (intrinsic) covariant derivative of the associated Phong-Rodrigues vector field //
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesCovariantDerivativeField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesCovariantDerivativeField< K , K+1 > , Point< double , K+1 > >;

		///////////////////////////////////////////////////////////////////////////////////////////////////
		// A class returning the intrinsic representation of the associated Phong-Rodrigues vector field //
		///////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesIntrinsicVectorField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesIntrinsicVectorField< K , K+1 > , Point< double , K+1 > >;

		/////////////////////////////////////////////////////////////////////////////////////
		// A class returning the divergence of the associated Phong-Rodrigues vector-field //
		/////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesDivergenceField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesDivergenceField< K , K+1 > , Point< double , K+1 > >;

		///////////////////////////////////////////////////////////////////////////////////
		// A class returning the bracket of the associated Phong-Rodrigues vector-fields //
		///////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesBracketField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesBracketField< K , K+1 > , Point< double , K+1 > , Point< double , K+1 > >;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// A class returning the covariant derivative of the one Phong-Rodrigues vector-field with respect to the other //
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		using PhongRodriguesCovariantDirectionalDerivativeField = EmbeddedPhongMeshField< K , SimplexProcessing::PhongRodriguesCovariantDirectionalDerivativeField< K , K+1 > , Point< double , K+1 > , Point< double , K+1 > >;

#include "SimplicialMeshProcessing.fields.inl"
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_FIELDS_INCLUDED
