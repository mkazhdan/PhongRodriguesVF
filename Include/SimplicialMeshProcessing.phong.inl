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

///////////////////////
// EmbeddedPhongMesh //
///////////////////////
template< unsigned int K >
EmbeddedPhongMesh< K >::EmbeddedPhongMesh( std::vector< Point< double , Dim > > &vertices , std::vector< Point< double , Dim > > &normals , std::vector< SimplexIndex< K > > &simplices )
	: _vertices(vertices) , _simplices(simplices) , _normals(normals)
{}

template< unsigned int K >
Simplex< double , EmbeddedPhongMesh< K >::Dim , K > EmbeddedPhongMesh< K >::simplexVertices( size_t s ) const
{
	Simplex< double , Dim , K > simplex;
	for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = _vertices[ _simplices[s][k] ];
	return simplex;
}

template< unsigned int K >
Simplex< double , EmbeddedPhongMesh< K >::Dim , K > EmbeddedPhongMesh< K >::simplexNormals( size_t s ) const
{
	Simplex< double , Dim , K > simplex;
	for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = _normals[ _simplices[s][k] ];
	return simplex;
}

template< unsigned int K >
auto EmbeddedPhongMesh< K >::metricTensorField( size_t sIdx ) const
{
	Point< double , Dim > v[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
	return [ g=SimplexProcessing::MetricTensorFromEmbedding< K , Dim >( v ) ]( SimplexProcessing::Position< K > ){ return g; };
}

template< unsigned int K >
auto EmbeddedPhongMesh< K >::measureScaleField( size_t sIdx ) const
{
	Point< double , Dim > v[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
	return [ scaleFactor=SimplexProcessing::ScaleFactorFromMetricTensor< K >( SimplexProcessing::MetricTensorFromEmbedding< K , Dim >( v ) ) ]( SimplexProcessing::Position< K > ){ return scaleFactor; };
}

template< unsigned int K >
auto EmbeddedPhongMesh< K >::inverseMetricTensorField( size_t sIdx ) const
{
	Point< double , Dim > v[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
	return [ gInv=SimplexProcessing::MetricTensorFromEmbedding< K , Dim >( v ).inverse() ]( SimplexProcessing::Position< K > ){ return gInv; };
}

template< unsigned int K >
auto EmbeddedPhongMesh< K >::normalField( size_t sIdx ) const
{
	Point< double , Dim > n[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) n[k] = _normals[ _simplices[sIdx][k] ];
	return SimplexProcessing::NormalizationField< K , Point< double , Dim > , SimplexProcessing::LinearInterpolant< K , Point< double , Dim > > >(n);
}

template< unsigned int K >
auto EmbeddedPhongMesh< K >::elementIndex( void ) const
{
	return [&]( size_t s , unsigned int e )
		{
			unsigned int k = e/(K+1) , _k = e%(K+1);
			return _simplices[s][k] * (K+1) + _k;
		};
}

template< unsigned int K >
auto EmbeddedPhongMesh< K >::elements( void ) const
{
	return SimplexProcessing::ArrayWrapper
	(
		[&]( size_t s )
		{
			Point< double , Dim > n[K+1];
			for( unsigned int k=0 ; k<=K ; k++ ) n[k] = _normals[ _simplices[s][k] ];
			return SimplexProcessing::PhongRodriguesSystem< K , Dim >::Elements( n ); 
		}
	);
}

template< unsigned int K >
template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , K > > , size_t > SampleFunctor >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::evaluationMatrix( size_t sampleNum , SampleFunctor && sampleFunctor ) const
{
	std::vector< Eigen::Triplet< double , size_t > > triplets( sampleNum * (K+1)*(K+1) * Dim );

	ThreadPool::ParallelFor
	(
		0 , sampleNum ,
		[&]( size_t s )
		{
			size_t idx = s * (K+1) * (K+1) * Dim;
			std::pair< size_t , Point< double , K > > sample = sampleFunctor( s );
			Point< double , Dim > n[K+1];
			for( unsigned int k=0 ; k<=K ; k++ ) n[k] = _normals[ _simplices[sample.first][k] ];
			typename SimplexProcessing::PhongRodriguesSystem< K , Dim >::Elements elements( n );

			for( unsigned int e=0 ; e<(K+1)*(K+1) ; e++ )
			{
				unsigned int k= e/(K+1) , _k = e%(K+1);
				Point< double , Dim > v = elements[e]( sample.second );
				for( unsigned int d=0 ; d<Dim ; d++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( s*Dim+d , _simplices[ sample.first ][k]*(K+1) + _k , v[d] );
			}
		}
	);

	Eigen::SparseMatrix< double > E( sampleNum * (K+1) , _vertices.size()*(K+1) );
	E.setFromTriplets( triplets.begin() , triplets.end() );
	return E;

}

template< unsigned int K >
#ifdef USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemVector > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
Eigen::VectorXd EmbeddedPhongMesh< K >::_systemVector( SystemVectorField && Sys , bool needsScaling , WeightField && WF ) const
{
	if constexpr( std::same_as< WeightField , UnitWeightField< K > > )
		return RiemannianMesh< K , EmbeddedPhongMesh< K > >::template _systemVector< QuadratureSamples , SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >( _vertices.size()*Dim , elementIndex() , std::forward< SystemVectorField >( Sys ) , needsScaling );
	else
	{
		auto WeightedSys = ScaledField< K >( std::forward< WeightField >( WF ) , std::forward< SystemVectorField >( Sys ) ); 
		return RiemannianMesh< K , EmbeddedPhongMesh< K > >::template _systemVector< QuadratureSamples , SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >( _vertices.size()*Dim , elementIndex() , WeightedSys , needsScaling );
	}
}

template< unsigned int K >
#ifdef USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemMatrix > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::_systemMatrix( SystemMatrixField && Sys , bool needsScaling , WeightField && WF ) const
{
	if constexpr( std::same_as< WeightField , UnitWeightField< K > > )
		return RiemannianMesh< K , EmbeddedPhongMesh< K > >::template _systemMatrix< QuadratureSamples , SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >( _vertices.size()*Dim , elementIndex() , std::forward< SystemMatrixField >( Sys ) , needsScaling );
	else
	{
		auto WeightedSys = ScaledField< K >( std::forward< WeightField >( WF ) , std::forward< SystemMatrixField >( Sys ) ); 
		return RiemannianMesh< K , EmbeddedPhongMesh< K > >::template _systemMatrix< QuadratureSamples , SimplexProcessing::PhongRodriguesSystem< K , Dim >::NumElements >( _vertices.size()*Dim , elementIndex() , WeightedSys , needsScaling );
	}
}

template< unsigned int K >
#ifdef USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemVector > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , Point< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemVectorField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
Eigen::VectorXd EmbeddedPhongMesh< K >::systemVector( SystemVectorField && Sys , WeightField && WF ) const
{
	return _systemVector< QuadratureSamples >( std::forward< SystemVectorField >( Sys ) , true , std::forward< WeightField >( WF ) );
}

template< unsigned int K >
#ifdef USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemMatrix > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
#else // !USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemMatrixField , HasMeshScaleFactorFunction< K > WeightField >
#endif // USING_GCC
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::systemMatrix( SystemMatrixField && Sys , WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( std::forward< SystemMatrixField >( Sys ) , true , std::forward< WeightField >( WF ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::Vector > VectorField , HasMeshScaleFactorFunction< K > WeightField >
Eigen::VectorXd EmbeddedPhongMesh< K >::massVector( VectorField && F , WeightField && WF ) const
{
	return _systemVector< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::PhongRodriguesSystem< K , Dim >::MassVector( simplexVertices(sIdx) , simplexNormals(sIdx) , F[sIdx] ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::massMatrix( WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::PhongRodriguesSystem< K , Dim >::MassMatrix( simplexVertices(sIdx) , simplexNormals(sIdx) ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::stiffnessMatrix( WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::PhongRodriguesSystem< K , Dim >::StiffnessMatrix( simplexVertices(sIdx) , simplexNormals(sIdx) ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , unsigned int Components , HasMeshScaleFactorFunction< K > WeightField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::stiffnessMatrix( WeightField && WF ) const
{
	return _systemMatrix< QuadratureSamples >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::PhongRodriguesSystem< K , Dim >::template ComponentStiffnessMatrix< Components >( simplexVertices(sIdx) , simplexNormals(sIdx) ); } ) , false , std::forward< WeightField >( WF ) );
}

template< unsigned int K >
#ifdef USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , typename EmbeddedPhongMesh< K >::SystemMatrix > SystemField >
#else // !USING_GCC
template< unsigned int QuadratureSamples , HasMeshFunction< K , SquareMatrix< double , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements > > SystemField >
#endif // USING_GCC
void EmbeddedPhongMesh< K >::setSystemMatrixEntries( EigenMatrixEntries & eme , SystemField && Sys ) const
{
	this->template _setSystemMatrixEntries< QuadratureSamples , SimplexProcessing::PhongRodriguesSystem< K , EmbeddedPhongMesh< K >::Dim >::NumElements >( eme , std::forward< SystemField >( Sys ) , true );
}

template< unsigned int K >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::J( void ) const requires( K==2 )
{
	auto NormalRotation = []( Point< double , Dim > n )
	{
		SquareMatrix< double , Dim > R;
		R(0,1) =  n[2] , R(1,0) = -n[2];
		R(0,2) = -n[1] , R(2,0) =  n[1];
		R(1,2) =  n[0] , R(2,1) = -n[0];

		for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) R(i,j) += n[i] * n[j];
		return R;
	};

	Eigen::SparseMatrix< double > J;
	std::vector< Eigen::Triplet< double > > triplets( _normals.size()*Dim*Dim );

	ThreadPool::ParallelFor
	(
		0 , _normals.size() ,
		[&]( size_t i )
		{
			SquareMatrix< double , Dim > R = NormalRotation( _normals[i] / Point< double , Dim >::Length( _normals[i] ) );
			for( unsigned int j=0 ; j<Dim ; j++ ) for( unsigned int k=0 ; k<Dim ; k++ )
				triplets[ i*Dim*Dim + j*Dim + k ] = Eigen::Triplet< double >( static_cast< int >( i*Dim+j ) , static_cast< int >( i*Dim+k ) , R(j,k) );
		}
	);

	J.resize( _normals.size()*Dim , _normals.size()*Dim );
	J.setFromTriplets( triplets.begin() , triplets.end() );
	return J;
}

template< unsigned int K >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::tangentProlongation( void ) const
{
	auto SetFrame = []( Point< double , Dim > n , Point< double , Dim > frame[K+1] )
	{
		frame[0] = n;

		for( unsigned int d=1 ; d<K ; d++ )
		{
			frame[d] = RandomSpherePoint< double , Dim >();
			while( true )
			{
				for( unsigned int dd=0 ; dd<d ; dd++ ) frame[d] -= Point< double , Dim >::Dot( frame[dd] , frame[d] ) * frame[dd];
				if( frame[d].squareNorm()>1e-10 )
				{
					frame[d] /= sqrt( frame[d].squareNorm() );
					break;
				}
			}
		}
		frame[K] = Point< double , Dim >::CrossProduct( frame );
	};

	Eigen::SparseMatrix< double > P;
	std::vector< Eigen::Triplet< double > > triplets( _normals.size()*Dim*K );

	ThreadPool::ParallelFor
	(
		0 , _normals.size() ,
		[&]( size_t i )
		{
			Point< double , Dim > n = _normals[i] / Point< double , Dim >::Length( _normals[i] );
			Point< double , Dim > t[K+1];
			SetFrame( n , t );

			for( unsigned int j=0 ; j<Dim ; j++ ) for( unsigned int k=0 ; k<K ; k++ )
				triplets[i*Dim*K+j*K+k] = Eigen::Triplet< double >( (int)i*Dim+j , (int)i*K+k , t[k+1][j] );
		}
	);

	P.resize( _normals.size()*Dim , _normals.size()*K );
	P.setFromTriplets( triplets.begin() , triplets.end() );
	return P;
}
