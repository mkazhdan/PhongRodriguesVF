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

#ifndef SIMPLICIAL_MESH_PROCESSING_PHONG_INCLUDED
#define SIMPLICIAL_MESH_PROCESSING_PHONG_INCLUDED

namespace MishaK
{
	namespace SimplicialMesh
	{
		enum struct CovariantComponent
		{
			Total ,
			Trace ,
			Traceless ,
			Symmetric ,
			AntiSymmetric ,
			Hodge
		};
		inline static const std::vector< std::string > CovariantComponentNames = { "total" , "trace" , "trace-less" , "symmetric" , "anti-symmetric" , "hodge" };

		template< unsigned int K >
		struct EmbeddedPhongMesh : public RiemannianMesh< K , EmbeddedPhongMesh< K > >
		{
			static const unsigned int Dim = K+1;
			using Scalar = double;
			using Vector = Point< double , Dim >;
			using EigenMatrixEntries = SystemIntegration::EigenMatrixEntries< SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements >;

			using RiemannianMesh< K , EmbeddedPhongMesh< K > >::       metricTensorField;
			using RiemannianMesh< K , EmbeddedPhongMesh< K > >::inverseMetricTensorField;
			using RiemannianMesh< K , EmbeddedPhongMesh< K > >::       measureScaleField;

			EmbeddedPhongMesh( std::vector< Point< double , Dim > > &vertices , std::vector< Point< double , Dim > > &normals , std::vector< SimplexIndex< K > > &simplices );

			operator EmbeddedMesh< K , K+1 >() { return EmbeddedMesh< K , K+1 >( _vertices , _simplices ); }
			operator const EmbeddedMesh< K , K+1 >() const { return EmbeddedMesh< K , K+1 >( _vertices , _simplices ); }

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
			Point< double , Dim > &normal( size_t n ){ return _normals[n]; }
			std::vector< Point< double , Dim > > &normals( void ){ return _normals; }

			const Point< double , Dim > &normal( size_t n ) const { return _normals[n]; }
			const std::vector< Point< double , Dim > > &normals( void ) const { return _normals; }

			Simplex< double , EmbeddedPhongMesh< K >::Dim , K > simplexVertices( size_t s ) const;
			Simplex< double , EmbeddedPhongMesh< K >::Dim , K > simplexNormals( size_t s ) const;

			auto        metricTensorField( size_t sIdx ) const;
			auto        measureScaleField( size_t sIdx ) const;
			auto inverseMetricTensorField( size_t sIdx ) const;
			auto              normalField( size_t sIdx ) const;

			auto scalarElementIndex( void ) const;
			auto vectorElementIndex( void ) const;

			auto scalarElements( void ) const;
			auto vectorElements( void ) const;

			template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , K > > , size_t > SampleFunctor >
			Eigen::SparseMatrix< double > scalarEvaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const;

			template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , K > > , size_t > SampleFunctor >
			Eigen::SparseMatrix< double > vectorEvaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const;


			template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::VectorXd scalarDual( ScalarOrDifferentialField && F , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField >
			Eigen::VectorXd scalarDual( ScalarOrDifferentialField && F ) const;

			template< unsigned int QuadratureSamples >
			Eigen::SparseMatrix< double > scalarMass( void ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > scalarMass( ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples >
			Eigen::SparseMatrix< double > scalarStiffness( void ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > scalarStiffness( ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , typename SystemField >
			Eigen::SparseMatrix< double > scalarSystem( SystemField && Sys ) const;

			template< unsigned int QuadratureSamples >
			SquareMatrix< double , K+1 > simplexScalarMass( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , K+1 > simplexScalarMass( size_t sIdx , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples >
			SquareMatrix< double , K+1 > simplexScalarStiffness( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , K+1 > simplexScalarStiffness( size_t sIdx , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , typename SystemField >
			SquareMatrix< double , K+1 > simplexScalarSystem( size_t sIdx , SystemField && Sys ) const;


			EigenMatrixEntries scalarEigenMatrixEntries( void ) const { return EigenMatrixEntries( _vertices.size() , _simplices.size() , scalarElementIndex() ); }

			template< unsigned int QuadratureSamples , typename SystemField >
			void setScalarSystemEntries( EigenMatrixEntries &eme , SystemField && Sys ) const;

			template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::VectorXd vectorDual( ScalarOrDifferentialField && F , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField >
			Eigen::VectorXd vectorDual( ScalarOrDifferentialField && F ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > vectorMass( ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples >
			Eigen::SparseMatrix< double > vectorMass( void ) const;

			template< unsigned int QuadratureSamples , typename SystemField >
			Eigen::SparseMatrix< double > vectorSystem( SystemField && Sys , bool needsScaling=true ) const;

			template< unsigned int QuadratureSamples >
			Eigen::SparseMatrix< double > vectorCovariantStiffness( void ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > vectorCovariantStiffness( ScaleFactorField && SF ) const;

			template< unsigned int QuadratureSamples , CovariantComponent CComponent >
			Eigen::SparseMatrix< double > vectorCovariantStiffness( void ) const;

			template< unsigned int QuadratureSamples , CovariantComponent CComponent , HasMeshScaleFactorFunction< K > ScaleFactorField >
			Eigen::SparseMatrix< double > vectorCovariantStiffness( ScaleFactorField && SF ) const;

			// Gives the squared-norm of the bracket with v
			template< unsigned int QuadratureSamples , typename TangentVectorField >
			Eigen::SparseMatrix< double > vectorBracketEnergy( TangentVectorField && V ) const;

			// Gives the squared-norm of the dot-product with v
			template< unsigned int QuadratureSamples , typename TangentVectorField >
			Eigen::SparseMatrix< double > vectorDotProductEnergy( TangentVectorField && V ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorMass( size_t sIdx , ScaleFactorField && S ) const;

			template< unsigned int QuadratureSamples >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorMass( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , typename SystemField >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorSystem( size_t sIdx , SystemField && Sys ) const;

			template< unsigned int QuadratureSamples >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorCovariantStiffness( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorCovariantStiffness( size_t sIdx , ScaleFactorField && SF ) const;

			template< unsigned int QuadratureSamples , CovariantComponent CComponent >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorCovariantStiffness( size_t sIdx ) const;

			template< unsigned int QuadratureSamples , CovariantComponent CComponent , HasMeshScaleFactorFunction< K > ScaleFactorField >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorCovariantStiffness( size_t sIdx , ScaleFactorField && SF ) const;

			// Gives the squared-norm of the bracket with v
			template< unsigned int QuadratureSamples , typename TangentVectorField >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorBracketEnergy( size_t sIdx , TangentVectorField && V ) const;

			// Gives the squared-norm of the dot-product with v
			template< unsigned int QuadratureSamples , typename TangentVectorField >
			SquareMatrix< double , (K+1)*(K+1) > simplexVectorDotProductEnergy( size_t sIdx , TangentVectorField && V ) const;


			EigenMatrixEntries vectorEigenMatrixEntries( void ) const { return this->template _eigenMatrixEntries< (K+1)*Dim >( _vertices.size()*Dim , vectorElementIndex() ); }

			template< unsigned int QuadratureSamples , typename SystemField >
			void setVectorSystemEntries( EigenMatrixEntries & eme , SystemField && Sys ) const;

			Eigen::SparseMatrix< double > vectorJ( void ) const;

			Eigen::SparseMatrix< double > vectorRotate( double radians ) const;

			Eigen::SparseMatrix< double > tangentProlongation( void ) const;

		protected:

			std::vector< Point< double , Dim > > &_vertices;
			std::vector< SimplexIndex< K > > &_simplices;

			std::vector< Point< double , Dim > > &_normals;
		};

#include "SimplicialMeshProcessing.phong.inl"
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_PHONG_INCLUDED
