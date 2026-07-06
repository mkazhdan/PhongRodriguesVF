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

#ifndef SIMPLICIAL_MESH_PROCESSING_EMBEDDED_INCLUDED
#define SIMPLICIAL_MESH_PROCESSING_EMBEDDED_INCLUDED

namespace MishaK
{
	namespace SimplicialMesh
	{
		template< unsigned int K , unsigned int Dim >
		struct EmbeddedMesh : public RiemannianMesh< K , EmbeddedMesh< K , Dim > >
		{
			using Scalar = double;
			using EigenMatrixEntries = SystemIntegration::EigenMatrixEntries< SimplexProcessing::LinearElements< K >::NumElements >;

			using RiemannianMesh< K , EmbeddedMesh< K , Dim > >::       metricTensorField;
			using RiemannianMesh< K , EmbeddedMesh< K , Dim > >::inverseMetricTensorField;
			using RiemannianMesh< K , EmbeddedMesh< K , Dim > >::       measureScaleField;

			EmbeddedMesh( std::vector< Point< double , Dim > > &vertices , std::vector< SimplexIndex< K > > &simplices ) : _vertices(vertices) , _simplices(simplices) {}

			Point< double , Dim > &vertex( size_t v ){ return _vertices[v]; }
			const Point< double , Dim > &vertex( size_t v ) const { return _vertices[v]; }

			std::vector< Point< double , Dim > > &vertices( void ){ return _vertices; }
			const std::vector< Point< double , Dim > > &vertices( void ) const { return _vertices; }

			size_t vertexNum( void ) const { return _vertices.size(); }

			SimplexIndex< K > &simplex( size_t s ){ return _simplices[s]; }
			const SimplexIndex< K > &simplex( size_t s ) const { return _simplices[s]; }

			std::vector< SimplexIndex< K > > &simplices( void ){ return _simplices; }
			const std::vector< SimplexIndex< K > > &simplices( void ) const { return _simplices; }

			size_t simplexNum( void ) const { return _simplices.size(); }

			Simplex< double , Dim , K > simplexVertices( size_t s ) const;

			template< typename T , SimplexProcessing::HasFunction< T , Point< double , Dim > > ImplicitFunction >
			auto implicitFunctionToScalarField( ImplicitFunction && F ) const;

			auto        metricTensorField( size_t sIdx ) const;
			auto        measureScaleField( size_t sIdx ) const;
			auto inverseMetricTensorField( size_t sIdx ) const;

			auto scalarElementIndex( void ) const;

			auto scalarElements( void ) const;

			template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , Dim > > , size_t > SampleFunctor >
			Eigen::SparseMatrix< double > evaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const;

			template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::VectorXd scalarDual( ScalarOrDifferentialField && F , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField >
			Eigen::VectorXd scalarDual( ScalarOrDifferentialField && F ) const
				requires SimplexProcessing::HasArrayOfSimplexFunctions< ScalarOrDifferentialField , K , Scalar > || SimplexProcessing::HasArrayOfSimplexFunctions< ScalarOrDifferentialField , K , SimplexProcessing::Differential< K , Scalar > >;

			template< unsigned int QuadratureSamples >
			Eigen::SparseMatrix< double > scalarMass( void ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > scalarMass( ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples >
			Eigen::SparseMatrix< double > scalarStiffness( void ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > scalarStiffness( ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , HasMeshTangentVectorFunction< K > TangentVectorField >
			Eigen::SparseMatrix< double > scalarDerivationSystem( TangentVectorField && VF ) const;

			template< unsigned int QuadratureSamples , typename SystemField >
			Eigen::SparseMatrix< double > scalarSystem( SystemField && Sys , bool needsScaling=true ) const;

			template< unsigned int QuadratureSamples >
			SquareMatrix< double , K+1 > simplexScalarMass( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , K+1 > simplexScalarMass( size_t sIdx , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples >
			SquareMatrix< double , K+1 > simplexScalarStiffness( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , K+1 > simplexScalarStiffness( size_t sIdx , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , HasMeshTangentVectorFunction< K > TangentVectorField >
			SquareMatrix< double , K+1 > simplexScalarDerivationSystem( size_t sIdx , TangentVectorField && VF ) const;

			template< unsigned int QuadratureSamples , typename SystemField >
			SquareMatrix< double , K+1 > simplexScalarSystem( size_t sIdx , SystemField && Sys ) const;


			EigenMatrixEntries scalarEigenMatrixEntries( void ) const { return this->template _eigenMatrixEntries< K+1 >( _vertices.size() , scalarElementIndex() ); }

			template< unsigned int QuadratureSamples , typename SystemField >
			void setScalarSystemEntries( EigenMatrixEntries &eme , SystemField && Sys ) const;

		protected:
			std::vector< Point< double , Dim > > &_vertices;
			std::vector< SimplexIndex< K > > &_simplices;
		};

#include "SimplicialMeshProcessing.embedded.inl"
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_EMBEDDED_INCLUDED
