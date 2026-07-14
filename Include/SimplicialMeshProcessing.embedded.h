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
#ifdef USING_GCC
			using Differential = SimplexProcessing::Differential< K , Scalar >;
#endif // USING_GCC
			using EigenMatrixEntries = SystemIntegration::EigenMatrixEntries< SimplexProcessing::ScalarSystem< K >::NumElements >;

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

			auto elementIndex( void ) const;
			auto elements( void ) const;

			template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , Dim > > , size_t > SampleFunctor >
			Eigen::SparseMatrix< double > evaluationMatrix( size_t sampleNum , SampleFunctor && sampleFunctor ) const;


			template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::VectorXd systemVector( SystemVectorField && Sys , WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > systemMatrix( SystemMatrixField && Sys , WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedMesh< K , Dim >::Scalar > ScalarField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::VectorXd massVector( ScalarField && F , WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshDifferentiableFunction< K , typename EmbeddedMesh< K , Dim >::Scalar > ScalarField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::VectorXd stiffnessVector( ScalarField && F , WeightField && WF=UnitWeightField< K >() ) const;

#ifdef USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedMesh< K , Dim >::Differential > DifferentialField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
#else // !USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , SimplexProcessing::Differential< K , typename EmbeddedMesh< K , Dim >::Scalar > > DifferentialField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
#endif // USING_GCC
			Eigen::VectorXd stiffnessVector( DifferentialField && F , WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > massMatrix( WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > stiffnessMatrix( WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshTangentVectorFunction< K > TangentVectorField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > derivationSystemMatrix( TangentVectorField && VF , WeightField && WF=UnitWeightField< K >() ) const;


			EigenMatrixEntries eigenMatrixEntries( void ) const { return this->template _eigenMatrixEntries< K+1 >( _vertices.size() , elementIndex() ); }

			template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemField >
			void setSystemMatrixEntries( EigenMatrixEntries & eme , SystemField && Sys ) const;

		protected:
			template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
			Eigen::VectorXd _systemVector( SystemVectorField && Sys , bool needsScaling , WeightField && WF ) const;

			template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
			Eigen::SparseMatrix< double > _systemMatrix( SystemMatrixField && Sys , bool needsScaling , WeightField && WF ) const;

			std::vector< Point< double , Dim > > &_vertices;
			std::vector< SimplexIndex< K > > &_simplices;
		};

#include "SimplicialMeshProcessing.embedded.inl"
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_EMBEDDED_INCLUDED
