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
auto EmbeddedMesh< K , Dim >::elementIndex( void ) const { return [&]( size_t s , unsigned int e ){ return _simplices[s][e]; }; }

template< unsigned int K , unsigned int Dim >
auto EmbeddedMesh< K , Dim >::elements( void ) const
{
	return SimplexProcessing::ArrayWrapper( [&]( size_t ){ return SimplexProcessing::ScalarSystem< K >::Elements(); } );
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
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::evaluationMatrix( size_t sampleNum , SampleFunctor && sampleFunctor ) const
{
	static typename SimplexProcessing::ScalarSystem< K >::Elements elements;

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
template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::VectorXd EmbeddedMesh< K , Dim >::_systemVector( SystemVectorField && Sys , bool needsScaling , WeightField && WF ) const
{
	if constexpr( std::same_as< WeightField , UnitWeightField< K > > )
		return RiemannianMesh< K , EmbeddedMesh< K , Dim > >::template _systemVector< QuadratureSamples , SimplexProcessing::ScalarSystem< K >::NumElements >( _vertices.size() , elementIndex() , std::forward< SystemVectorField >( Sys ) , needsScaling );
	else
	{
		auto WeightedSys = ScaledField< K >( std::forward< WeightField >( WF ) , std::forward< SystemVectorField >( Sys ) ); 
		return RiemannianMesh< K , EmbeddedMesh< K , Dim > >::template _systemVector< QuadratureSamples , SimplexProcessing::ScalarSystem< K >::NumElements >( _vertices.size() , elementIndex() , WeightedSys , needsScaling );
	}
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::_systemMatrix( SystemMatrixField && Sys , bool needsScaling , WeightField && WF ) const
{
	if constexpr( std::same_as< WeightField , UnitWeightField< K > > )
		return RiemannianMesh< K , EmbeddedMesh< K , Dim > >::template _systemMatrix< QuadratureSamples , SimplexProcessing::ScalarSystem< K >::NumElements >( _vertices.size() , elementIndex() , std::forward< SystemMatrixField >( Sys ) , needsScaling );
	else
	{
		auto WeightedSys = ScaledField< K >( std::forward< WeightField >( WF ) , std::forward< SystemMatrixField >( Sys ) ); 
		return RiemannianMesh< K , EmbeddedMesh< K , Dim > >::template _systemMatrix< QuadratureSamples , SimplexProcessing::ScalarSystem< K >::NumElements >( _vertices.size() , elementIndex() , WeightedSys , needsScaling );
	}
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::VectorXd EmbeddedMesh< K , Dim >::systemVector( SystemVectorField && Sys , WeightField && WF ) const
{
	return _systemVector< QuadratureSamples >( std::forward< SystemVectorField >( Sys ) , true , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::systemMatrix( SystemMatrixField && Sys , WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( std::forward< SystemMatrixField >( Sys ) , true , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedMesh< K , Dim >::Scalar > ScalarField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::VectorXd EmbeddedMesh< K , Dim >::massVector( ScalarField && F , WeightField && WF ) const
{
	return _systemVector< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::ScalarSystem< K >::MassVector( simplexVertices( sIdx ) , F[sIdx] ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshDifferentiableFunction< K , typename EmbeddedMesh< K , Dim >::Scalar > ScalarField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::VectorXd EmbeddedMesh< K , Dim >::stiffnessVector( ScalarField && F , WeightField && WF ) const
{
	return _systemVector< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::ScalarSystem< K >::StiffnessVector( simplexVertices( sIdx ) , F[sIdx] ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
#ifdef USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedMesh< K , Dim >::Differential > DifferentialField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , SimplexProcessing::Differential< K , typename EmbeddedMesh< K , Dim >::Scalar > > DifferentialField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
Eigen::VectorXd EmbeddedMesh< K , Dim >::stiffnessVector( DifferentialField && F , WeightField && WF ) const
{
	return _systemVector< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::ScalarSystem< K >::StiffnessVector( simplexVertices( sIdx ) , F[sIdx] ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::massMatrix( WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::ScalarSystem< K >::MassMatrix( simplexVertices( sIdx ) ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::stiffnessMatrix( WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::ScalarSystem< K >::StiffnessMatrix( simplexVertices( sIdx ) ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshTangentVectorFunction< K > TangentVectorField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedMesh< K , Dim >::derivationSystemMatrix( TangentVectorField && VF , WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::ScalarSystem< K >::DerivationMatrix( simplexVertices( sIdx ) , VF[sIdx] ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K , unsigned int Dim >
template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::ScalarSystem< K >::NumElements > > SystemField >
void EmbeddedMesh< K , Dim >::setSystemMatrixEntries( EigenMatrixEntries &eme , SystemField && Sys ) const
{
	this->template _setSystemMatrixEntries< QuadratureSamples , SimplexProcessing::ScalarSystem< K >::NumElements >( eme , std::forward< SystemField >( Sys ) , true );
}
