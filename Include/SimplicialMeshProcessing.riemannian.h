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

#ifndef SIMPLICIAL_MESH_PROCESSING_RIEMANNIAN_INCLUDED
#define SIMPLICIAL_MESH_PROCESSING_RIEMANNIAN_INCLUDED

namespace MishaK
{
	namespace SimplicialMesh
	{
		// A general Reimannian mesh wrapper, use to enclose a mesh type that:
		// -- Defines a metric tensor function
		template< unsigned int K , typename MeshType >
		struct RiemannianMesh
		{
			// <---- Forwarding functionality from MeshType
			      SimplexIndex< K > &simplex( size_t s )       { return static_cast<       MeshType * >( this )->simplex(s); }
			const SimplexIndex< K > &simplex( size_t s ) const { return static_cast< const MeshType * >( this )->simplex(s); }

			      std::vector< SimplexIndex< K > > &simplices( void )       { return static_cast<       MeshType * >( this )->simplices(); }
			const std::vector< SimplexIndex< K > > &simplices( void ) const { return static_cast< const MeshType * >( this )->simplices(); }

			size_t simplexNum( void ) const { return static_cast< const MeshType * >( this )->simplexNum(); }

			auto metricTensorField( size_t sIdx ) const { return static_cast< const MeshType * >( this )->metricTensorField( sIdx ); }
			auto inverseMetricTensorField( size_t sIdx ) const 
			{
				if constexpr( _HasInverseMetricTensorField< MeshType >::value ) return static_cast< const MeshType * >( this )->inverseMetricTensorField( sIdx );
				else return [g=this->metricTensorField(sIdx)]( SimplexProcessing::Position< K > p ){ return g(p).inverse(); };
			}
			auto measureScaleField( size_t sIdx ) const 
			{
				if constexpr( _HasMeasureScaleField< MeshType >::value ) return static_cast< const MeshType * >( this )->measureScaleField( sIdx );
				else return [g=this->metricTensorField(sIdx)]( SimplexProcessing::Position< K > p ){ return sqrt( g(p).determinant() ); };
			}
			//  Forwarding functionality from MeshType ---->

			auto        metricTensorField( void ) const { return SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return this->       metricTensorField(sIdx); } ); }
			auto inverseMetricTensorField( void ) const { return SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return this->inverseMetricTensorField(sIdx); } ); }
			auto        measureScaleField( void ) const { return SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return this->       measureScaleField(sIdx); } ); }

			template< unsigned int QuadratureSamples >
			double measure( void ) const;

			template< unsigned int QuadratureSamples , typename T , HasMeshFunction< K , T > ScalarField >
			T integral( ScalarField && F ) const;

		protected:
			template< typename T , typename=void >
			struct _HasInverseMetricTensorField : std::false_type {};

			template< typename T >
			struct _HasInverseMetricTensorField< T , std::void_t< decltype( std::declval< T >().inverseMetricTensorField( std::declval< size_t >() ) ) > > : std::true_type {};

			template< typename T , typename=void >
			struct _HasMeasureScaleField : std::false_type {};

			template< typename T >
			struct _HasMeasureScaleField< T , std::void_t< decltype( std::declval< T >().measureScaleField( std::declval< size_t >() ) ) > > : std::true_type {};


			template< typename T , HasMeshScaleFactorFunction< K > ScaleFactorField /* = SimplexProcessing::Samples< SimplexProcessing::ScaleFactor > */ >
			auto _scaledIdentityField( ScaleFactorField && S ) const;

			template< typename T , HasMeshScaleFactorFunction< K > ScaleFactorField >
			auto _scaledInverseMetricTensorField( ScaleFactorField && S ) const;

			template< typename T , HasDifferentiableElements< K , T > Elements >
			static auto _DifferentialElements( Elements && E );

			template< typename T , HasElementsAndDifferentiableElements< K , T > Elements >
			static auto _ValueAndDifferentialElements( Elements && E );


			// Functionality for defining system vertices/matrices 
			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshFunction< K , T > ValueField , HasElements< K , T > Elements , HasElementIndexFunctor ElementIndex >
			Eigen::VectorXd _dual( size_t fNum , Elements && E , ValueField && F , ElementIndex && Idx ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshDifferentialFunction< K , T > DifferentialField , HasDifferentiableElements< K , T > Elements , HasElementIndexFunctor ElementIndex >
			Eigen::VectorXd _dual( size_t fNum , Elements && E , DifferentialField && F , ElementIndex && Idx ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasElements< K , T > Elements , HasElementIndexFunctor ElementIndex >
			Eigen::SparseMatrix< double > _mass( size_t fNum , Elements && E , ElementIndex && Idx ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasDifferentiableElements< K , T > Elements , HasElementIndexFunctor ElementIndex >
			Eigen::SparseMatrix< double > _stiffness( size_t fNum , Elements && E , ElementIndex && Idx ) const;


			// Functionality for defining weighted system vertices/matrices 
			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshFunction< K , T > ValueField , HasElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshScaleFactorFunction< K > WeightField >
			Eigen::VectorXd _weightedDual( size_t fNum , Elements && E , ValueField && F , ElementIndex && Idx , WeightField && WF ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshDifferentialFunction< K , T > DifferentialField , HasDifferentiableElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshScaleFactorFunction< K > WeightField >
			Eigen::VectorXd _weightedDual( size_t fNum , Elements && E , DifferentialField && F , ElementIndex && Idx , WeightField && WF ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshScaleFactorFunction< K > WeightField >
			Eigen::SparseMatrix< double > _weightedMass( size_t fNum , Elements && E , ElementIndex && Idx , WeightField && WF ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasDifferentiableElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshScaleFactorFunction< K > WeightField >
			Eigen::SparseMatrix< double > _weightedStiffness( size_t fNum , Elements && E , ElementIndex && Idx , WeightField && WF ) const;


			// Functionality for defining system matrices 
			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshSystemLinearMapOrBilinearFormFunction< K , NumElementsPerSimplex , T > SystemField >
			Eigen::SparseMatrix< double > _system( size_t fNum , Elements && E , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasDifferentiableElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
			Eigen::SparseMatrix< double > _system( size_t fNum , Elements && E , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasElementsAndDifferentiableElements< K , T > Elements , HasElementIndexFunctor ElementIndex , HasMeshSystemAndDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
			Eigen::SparseMatrix< double > _system( size_t fNum , Elements && E , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , HasElementIndexFunctor ElementIndex , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > SystemField >
			Eigen::SparseMatrix< double > _system( size_t fNum ,                 ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const;


			// Functionality for setting/resetting matrix entries
			template< unsigned int NumElementsPerSimplex , HasElementIndexFunctor ElementIndex >
			SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > _eigenMatrixEntries( size_t fNum , ElementIndex && Idx ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasElements< K , T > Elements , HasMeshSystemLinearMapOrBilinearFormFunction< K , NumElementsPerSimplex , T > SystemField >
			void _setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , Elements && E , SystemField && Sys ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasDifferentiableElements< K , T > Elements , HasMeshDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
			void _setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , Elements && E , SystemField && Sys ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasElementsAndDifferentiableElements< K , T > Elements , HasMeshSystemAndDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
			void _setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , Elements && E , SystemField && Sys ) const;

			template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > SystemField >
			void _setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme ,                 SystemField && Sys ) const;


			template< unsigned int NumElementsPerSimplex , typename T , HasMeshTangentVectorFunction< K > TangentVectorField >
			static auto _DerivationSystemField( TangentVectorField && VF );
		};
#include "SimplicialMeshProcessing.riemannian.inl"
	}
}

#endif // SIMPLICIAL_MESH_PROCESSING_RIEMANNIAN_INCLUDED
