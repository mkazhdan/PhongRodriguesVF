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
		template< unsigned int K >
		struct EmbeddedPhongMesh : public RiemannianMesh< K , EmbeddedPhongMesh< K > >
		{
			static const unsigned int Dim = K+1;
			using Scalar = double;
			using Vector = Point< double , Dim >;
			using EigenMatrixEntries = SystemIntegration::EigenMatrixEntries< SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >;
#ifdef USING_GCC
			using SystemVector = Point< double , SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >;
			using SystemMatrix = SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >;
#endif // USING_GCC

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

			auto elementIndex( void ) const;
			auto elements( void ) const;

			template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , K > > , size_t > SampleFunctor >
			Eigen::SparseMatrix< double > evaluationMatrix( size_t sampleNum , SampleFunctor && sampleFunctor ) const;


#ifdef USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemVector > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
#else // !USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
#endif // USING_GCC
			Eigen::VectorXd systemVector( SystemVectorField && Sys , WeightField && WF=UnitWeightField< K >() ) const;

#ifdef USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemMatrix > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
#else // !USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
#endif // USING_GCC
			Eigen::SparseMatrix< double > systemMatrix( SystemMatrixField && Sys , WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::Vector > VectorField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::VectorXd massVector( VectorField && F , WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > massMatrix( WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > stiffnessMatrix( WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , unsigned int Components , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > stiffnessMatrix( WeightField && WF=UnitWeightField< K >() ) const;

			template< unsigned int QuadratureSamples , HasMeshDifferentiableFunction< K , typename EmbeddedPhongMesh< K >::Vector > VectorField , HasMeshScaleFactorFunction< K > WeightField=UnitWeightField< K > >
			Eigen::SparseMatrix< double > lieBracketMassMatrix( VectorField && VF , WeightField && WF=UnitWeightField< K >() ) const;


			EigenMatrixEntries eigenMatrixEntries( void ) const { return this->template _eigenMatrixEntries< (K+1)*Dim >( _vertices.size()*Dim , elementIndex() ); }

#ifdef USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemMatrix > SystemField >
#else // !USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemField >
#endif // USING_GCC
			void setSystemMatrixEntries( EigenMatrixEntries & eme , SystemField && Sys ) const;

			Eigen::SparseMatrix< double > J( void ) const requires( K==2 );

			Eigen::SparseMatrix< double > tangentProlongation( void ) const;

		protected:
#ifdef USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemVector > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
			Eigen::VectorXd _systemVector( SystemVectorField && Sys , bool needsScaling , WeightField && WF ) const;

#ifdef USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemMatrix > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
			template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
			Eigen::SparseMatrix< double > _systemMatrix( SystemMatrixField && Sys , bool needsScaling , WeightField && WF ) const;

			std::vector< Point< double , Dim > > &_vertices;
			std::vector< SimplexIndex< K > > &_simplices;

			std::vector< Point< double , Dim > > &_normals;
		};

#include "SimplicialMeshProcessing.phong.inl"
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_PHONG_INCLUDED
