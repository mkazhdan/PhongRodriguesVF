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

#ifndef SIMPLICIAL_MESH_PROCESSING_INTEGRATION_INCLUDED
#define SIMPLICIAL_MESH_PROCESSING_INTEGRATION_INCLUDED

namespace MishaK
{
	namespace SimplicialMesh
	{
		namespace SystemIntegration
		{
			template< unsigned int K , unsigned int NumF , typename T , typename Function /* = Field< K , T > */ , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
			auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B );

			template< unsigned int K , unsigned int NumF , typename T , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
			auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B );

			template< unsigned int K , unsigned int NumF , typename T , typename ScaleFactors /* = Samples< Field< K , ScaleFactor > > */ , typename Function /* = Field< K , T > */ , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
			auto ScaledSystemVectorField( ScaleFactors && SF , Function && F , TestFunctions && Fs , BilinearForms && B );

			template< unsigned int K , unsigned int NumF , typename T , typename ScaleFactors /* = Samples< Field< K , ScaleFactor > > */ , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
			auto ScaledSystemMatrixField( ScaleFactors && SF , TestFunctions && Fs , BilinearForms && B );

			template< unsigned int K , unsigned int NumF , typename T , typename ScaleFactors /* = Samples< Field< K , ScaleFactor > > */ , typename Field /* = Samples< K , T > */ >
			auto ScaledField( ScaleFactors && SF , Field && F );

			template< unsigned int NumElementsPerSimplex >
			struct EigenMatrixEntries
			{
				EigenMatrixEntries( void ) : _M( new Eigen::SparseMatrix< double >() ) {}
				~EigenMatrixEntries( void ){ delete _M; }

				EigenMatrixEntries( EigenMatrixEntries && e ){ std::swap(_M,e._M) , std::swap( _matrixEntries , e._matrixEntries ); }
				EigenMatrixEntries & operator = ( EigenMatrixEntries && e ){ std::swap(_M,e._M) , std::swap( _matrixEntries , e._matrixEntries ) ; return *this; }

				template< HasIndexFunctor Index >
				EigenMatrixEntries( size_t numF , size_t numS , Index && Idx );

				void clear( void );
				size_t size( void ) const { return _matrixEntries.size(); }
				double &matrixEntry( size_t i ){ return *_matrixEntries[i]; }

				const Eigen::SparseMatrix< double > &operator()( void ) const { return *_M; }
			protected:
				Eigen::SparseMatrix< double > *_M;
				std::vector< double * > _matrixEntries;

				template< unsigned int _NumElementsPerSimplex , typename _MassFunctor /* Samples< SquareMatrix< double , NumElementsPerSimplex > > */ >
				friend void _SetMassEntries( EigenMatrixEntries< _NumElementsPerSimplex > &eme , unsigned int numS , _MassFunctor && Mass );
			};

			template< unsigned int K , unsigned int QuadratureSamples >
			struct MCIntegrator
			{
				// The integral of a function over the simplex
				template< typename T , HasMeshFunction< K , T > Function >
				static T Integral( size_t numS , Function && F );

				// The integral of one function against a set of test functions over the simplex
				template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , HasMeshFunction< K , Point< double , NumElementsPerSimplex > > VectorFunctor >
				static Eigen::VectorXd Vector( size_t numF , size_t numS , Index && Idx , VectorFunctor && VF );

				// The integral of a set of test functions against each other over the simplex
				template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
				static Eigen::SparseMatrix< double > Matrix( size_t numF , size_t numS , Index && Idx , MatrixFunctor && MF );

				// The integral of a set of test functions against each other over the simplex
				template< unsigned int NumElementsPerSimplex , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
				static void SetMatrixEntries( EigenMatrixEntries< NumElementsPerSimplex > &eme , size_t numS , MatrixFunctor && MF );
			};

			//////////////////////////
			// FEM assembly helpers //
			//////////////////////////
			struct SystemAssembler
			{
				template< typename T , SimplexProcessing::HasArray< T > Integrand >
				static T Integral( size_t numS , Integrand && integrand );

				template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , SimplexProcessing::HasArray< Point< double , NumElementsPerSimplex > > VectorFunctor >
				static Eigen::VectorXd Vector( size_t numF , size_t numS , Index && Idx , VectorFunctor && VF );

				template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , SimplexProcessing::HasArray< SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
				static Eigen::SparseMatrix< double > Matrix( size_t numF , size_t numS , Index && Idx , MatrixFunctor && MF );

				template< unsigned int NumElementsPerSimplex , SimplexProcessing::HasArray< SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
				static void SetMatrixEntries( EigenMatrixEntries< NumElementsPerSimplex > &eme , size_t numS , MatrixFunctor && MF );
			};

#include "SimplicialMeshProcessing.integration.inl"
		}
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_INTEGRATION_INCLUDED
