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

//////////////////
// EmbeddedMesh //
//////////////////

template< unsigned int K , unsigned int Dim >
Simplex< double , Dim , K > EmbeddedMesh< K , Dim >::simplexVertices( size_t s ) const
{
	Simplex< double , Dim , K > simplex;
	for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = _vertices[ _simplices[s][k] ];
	return simplex;
}

template< unsigned int K , unsigned int Dim >
auto EmbeddedMesh< K , Dim >::metricTensorField( size_t sIdx ) const
{
	Point< double , Dim > v[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
	return [ g=SimplexProcessing::MetricTensorFromEmbedding< K , Dim >( v ) ]( SimplexProcessing::Position< K > ){ return g; };
}

template< unsigned int K , unsigned int Dim >
auto EmbeddedMesh< K , Dim >::measureScaleField( size_t sIdx ) const
{
	Point< double , Dim > v[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
	return [ scaleFactor=SimplexProcessing::ScaleFactorFromMetricTensor< K >( SimplexProcessing::MetricTensorFromEmbedding< K , Dim >( v ) ) ]( SimplexProcessing::Position< K > ){ return scaleFactor; };
}

template< unsigned int K , unsigned int Dim >
auto EmbeddedMesh< K , Dim >::inverseMetricTensorField( size_t sIdx ) const
{
	Point< double , Dim > v[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
	return [ gInv=SimplexProcessing::MetricTensorFromEmbedding< K , Dim >( v ).inverse() ]( SimplexProcessing::Position< K > ){ return gInv; };
}

template< unsigned int K , unsigned int Dim >
auto EmbeddedMesh< K , Dim >::scalarElementIndex( void ) const { return [&]( size_t s , unsigned int e ){ return _simplices[s][e]; }; }

template< unsigned int K , unsigned int Dim >
auto EmbeddedMesh< K , Dim >::scalarElements( void ) const
{
	return SimplexProcessing::ArrayWrapper( [&]( size_t ){ return SimplexProcessing::LinearElements< K >(); } );
}

template< unsigned int K , unsigned int Dim >
template< typename T , SimplexProcessing::HasFunction< T , Point< double , Dim > > ImplicitFunction >
auto EmbeddedMesh< K , Dim >::implicitFunctionToScalarField( ImplicitFunction && F ) const
{
	auto _F = [&]( size_t sIdx )
	{
		Simplex< double , Dim , K > s;
		for( unsigned int k=0 ; k<=K ; k++ ) s[k] = _vertices[ _simplices[sIdx][k] ];
		return [F,s]( SimplexProcessing::Position< K > p ){ return F( s(p) ); };
	};
	return SimplexProcessing::ArrayWrapper( _F );
}

template< unsigned int K , unsigned int Dim >
template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , Dim > > , size_t > SampleFunctor >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::evaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const
{
	static SimplexProcessing::LinearElements< K > elements;

	std::vector< Eigen::Triplet< double , size_t > > triplets( sampleNum * (K+1) );

	ThreadPool::ParallelFor
	(
		0 , sampleNum ,
		[&]( size_t s )
		{
			size_t idx = s * (K+1);
			std::pair< size_t , Point< double , K > > sample = sampleFunctor( s );

			for( unsigned int k=0 ; k<=K ; k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( s , _simplices[ sample.first ][k] , elements[k]( sample.second ) );
		}
	);

	Eigen::SparseMatrix< double > E( sampleNum , _vertices.size() );
	E.setFromTriplets( triplets.begin() , triplets.end() );
	return E;
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::VectorXd EmbeddedMesh< K , Dim >::scalarDual( ScalarOrDifferentialField && F , ScaleFactorField && S ) const
{
	return this->template _dual< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , std::forward< ScalarOrDifferentialField >( F ) , scalarElementIndex() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField >
Eigen::VectorXd EmbeddedMesh< K , Dim >::scalarDual( ScalarOrDifferentialField && F ) const
	requires SimplexProcessing::HasArrayOfSimplexFunctions< ScalarOrDifferentialField , K , Scalar > || SimplexProcessing::HasArrayOfSimplexFunctions< ScalarOrDifferentialField , K , SimplexProcessing::Differential< K , Scalar > >
{
	return this->template _dual< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , std::forward< ScalarOrDifferentialField >( F ) , scalarElementIndex() );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::scalarMass( void ) const
{
	return this->template _mass< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , scalarElementIndex() );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::scalarMass( ScaleFactorField && S ) const
{
	return this->template _mass< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , scalarElementIndex() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::scalarStiffness( void ) const
{
	return this->template _stiffness< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , scalarElementIndex() );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::scalarStiffness( ScaleFactorField && S ) const
{
	return this->template _stiffness< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , scalarElementIndex() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshTangentVectorFunction< K > TangentVectorField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::scalarDerivationSystem( TangentVectorField && VF ) const
{
	return this->template _system< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , scalarElementIndex() , RiemannianMesh< K , EmbeddedMesh< K , Dim > >::template _DerivationSystemField< K+1 , Scalar >( VF ) , true );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , typename SystemField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::scalarSystem( SystemField && Sys , bool needsScaling ) const
{
	return this->template _system< QuadratureSamples , K+1 , Scalar >( _vertices.size() , scalarElements() , scalarElementIndex() , std::forward< SystemField >( Sys ) , needsScaling );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples >
SquareMatrix< double , K+1 > EmbeddedMesh< K , Dim >::simplexScalarMass( size_t sIdx ) const
{
	return this->template _simplexMmss< QuadratureSamples , K+1 , Scalar >( sIdx , scalarElements() );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
SquareMatrix< double , K+1 > EmbeddedMesh< K , Dim >::simplexScalarMass( size_t sIdx , ScaleFactorField && S ) const
{
	return this->template _simplexMmss< QuadratureSamples , K+1 , Scalar >( sIdx , scalarElements() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples >
SquareMatrix< double , K+1 > EmbeddedMesh< K , Dim >::simplexScalarStiffness( size_t sIdx ) const
{
	return this->template _simplexStiffness< QuadratureSamples , K+1 , Scalar >( sIdx , scalarElements() );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
SquareMatrix< double , K+1 > EmbeddedMesh< K , Dim >::simplexScalarStiffness( size_t sIdx , ScaleFactorField && S ) const
{
	return this->template _simplexStiffness< QuadratureSamples , K+1 , Scalar >( sIdx , scalarElements() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshTangentVectorFunction< K > TangentVectorField >
SquareMatrix< double , K+1 > EmbeddedMesh< K , Dim >::simplexScalarDerivationSystem( size_t sIdx , TangentVectorField && VF ) const
{
	return this->template _simplexSystem< QuadratureSamples , K+1 , Scalar >( sIdx , scalarElements() , RiemannianMesh< K , EmbeddedMesh< K , Dim > >::template _DerivationSystemField< K+1 , Scalar >( VF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , typename SystemField >
SquareMatrix< double , K+1 > EmbeddedMesh< K , Dim >::simplexScalarSystem( size_t sIdx , SystemField && Sys ) const
{
	return this->template _simplexSystem< QuadratureSamples , K+1 , Scalar >( sIdx , scalarElements() , std::forward< SystemField >( Sys ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , typename SystemField >
void EmbeddedMesh< K , Dim >::setScalarSystemEntries( EigenMatrixEntries &eme , SystemField && Sys ) const
{
	this->template _setSystemEntries< QuadratureSamples , K+1 , Scalar >( eme , scalarElements() , std::forward< SystemField >( Sys ) );
}
