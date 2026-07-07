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
			return SimplexProcessing::PhongRodriguesVectorElements< K , Dim >( n ); 
		}
	);
}

template< unsigned int K >
template< SimplexProcessing::HasFunction< std::pair< size_t , Point< double , K > > , size_t > SampleFunctor >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::evaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const
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
			SimplexProcessing::PhongRodriguesVectorElements< K , Dim > elements( n );

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
template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::VectorXd EmbeddedPhongMesh< K >::dual( ScalarOrDifferentialField && F , ScaleFactorField && S ) const
{
	return this->template _dual< QuadratureSamples , (K+1)*Dim , Vector >( _vertices.size()*Dim , elements() , std::forward< ScalarOrDifferentialField >( F ) , elementIndex() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , typename ScalarOrDifferentialField >
Eigen::VectorXd EmbeddedPhongMesh< K >::dual( ScalarOrDifferentialField && F ) const
{
	return this->template _dual< QuadratureSamples , (K+1)*Dim , Vector >( _vertices.size()*Dim , elements() , std::forward< ScalarOrDifferentialField >( F ) , elementIndex() );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::mass( ScaleFactorField && S ) const
{
	return this->template _mass< QuadratureSamples , (K+1)*Dim , Vector >( _vertices.size()*Dim , elements() , elementIndex() , std::forward< ScaleFactorField >( S ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::mass( void ) const
{
	return this->template _mass< QuadratureSamples , (K+1)*Dim , Vector >( _vertices.size()*Dim , elements() , elementIndex() );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , typename SystemField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::system( SystemField && Sys , bool needsScaling ) const
{
	return this->template _system< QuadratureSamples , (K+1)*Dim , Vector >( _vertices.size()*Dim , elements() , elementIndex() , std::forward< SystemField >( Sys ) , needsScaling );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , typename SystemField >
void EmbeddedPhongMesh< K >::setSystemEntries( EigenMatrixEntries &eme , SystemField && Sys ) const
{
	this->template _setSystemEntries< QuadratureSamples , (K+1)*Dim , Vector >( eme , elements() , std::forward< SystemField >( Sys ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::stiffness( void ) const
{
	// Represent the covariant derivative as a map from TM -> T^*M
	static const unsigned int NumE = SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements;

	auto Sys = [&]( size_t sIdx )
		{
			return [tN=this->normalField(sIdx),tG=this->metricTensorField(sIdx),tS=this->measureScaleField(sIdx)]( SimplexProcessing::Position< K > p )
				{
					return [tn=tN(p),g=tG(p),s=tS(p)]( const SimplexProcessing::Differential< K , Vector > D[] )
						{
							SquareMatrix< double , K > _gAdj = g.adjugate() / s;
							Matrix< double , K , Dim > cov[ NumE ];
							for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int k=0 ; k<K ; k++ )
							{
								Point< double , Dim > _D = D[n][k] - tn * Point< double , Dim >::Dot( D[n][k] , tn );
								for( unsigned int d=0 ; d<Dim ; d++ ) cov[n](k,d) = _D[d];
							}

							SquareMatrix< double , NumE > mass;
							// The inner product of two maps M,N: TM -> R^d is
							//		<M,N> = tr( ( gInv * N.transpose() * M );
							for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int m=0 ; m<NumE ; m++ )
								mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj * cov[n].transpose() * cov[m] );
							return mass;
						};
				};
		};
	return this->template system< QuadratureSamples >( SimplexProcessing::ArrayWrapper( Sys ) , false );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::stiffness( ScaleFactorField && SF ) const
{
	// Represent the covariant derivative as a map from TM -> T^*M
	static const unsigned int NumE = SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements;

	auto Sys = [&]( size_t sIdx )
		{
			return [tN=this->normalField(sIdx),tG=this->metricTensorField(sIdx),tS=this->measureScaleField(sIdx),tSF=SF[sIdx]]( SimplexProcessing::Position< K > p )
				{
					return [tn=tN(p),g=tG(p),s=tS(p),sf=tSF(p)]( const SimplexProcessing::Differential< K , Vector > D[] )
						{
							SquareMatrix< double , K > _gAdj = g.adjugate() / s;
							Matrix< double , K , Dim > cov[ NumE ];
							for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int k=0 ; k<K ; k++ )
							{
								Point< double , Dim > _D = D[n][k] - tn * Point< double , Dim >::Dot( D[n][k] , tn );
								for( unsigned int d=0 ; d<Dim ; d++ ) cov[n](k,d) = _D[d];
							}

							SquareMatrix< double , NumE > mass;
							// The inner product of two maps M,N: TM -> R^d is
							//		<M,N> = tr( ( gInv * N.transpose() * M );
							for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int m=0 ; m<NumE ; m++ )
								mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj * cov[n].transpose() * cov[m] );
							return mass * sf;
						};
				};
		};
	return this->template system< QuadratureSamples >( SimplexProcessing::ArrayWrapper( Sys ) , false );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , CovariantComponent CComponent >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::stiffness( void ) const
{
	// Represent the covariant derivative as a map from TM -> T^*M
	static const unsigned int NumE = SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements;

	auto Sys = [&]( size_t sIdx )
	{
		return [tI2E=SimplexProcessing::PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >( simplexVertices(sIdx) , simplexNormals(sIdx) ),tG=this->metricTensorField(sIdx),tS=this->measureScaleField(sIdx)]( SimplexProcessing::Position< K > p )
		{
				return [i2e=tI2E(p),g=tG(p),s=tS(p)]( const SimplexProcessing::Differential< K , Vector > d[] )
			{
				SquareMatrix< double , K > _gAdj = g.adjugate() / s , _g = g / s;
				SquareMatrix< double , K > cov[ NumE ] , sym[ NumE ] , aSym[ NumE ] , trace[ NumE ];
				Point< double , K > e;
				for( unsigned int n=0 ; n<NumE ; n++ )
				{
					for( unsigned int k=0 ; k<K ; k++ )
					{
						e[k] = 1.;
						Point< double , Dim > _d = i2e * e;
						for( unsigned int _k=0 ; _k<K ; _k++ ) cov[n](k,_k) = Point< double , Dim >::Dot( _d , d[n][_k] );
						e[k] = 0.;
					}

					if constexpr( CComponent==CovariantComponent::AntiSymmetric || CComponent==CovariantComponent::Hodge     )  aSym[n] = ( cov[n] - cov[n].transpose() ) / 2.;
					if constexpr( CComponent==CovariantComponent::Symmetric     || CComponent==CovariantComponent::Traceless )   sym[n] = ( cov[n] + cov[n].transpose() ) / 2.;
					if constexpr( CComponent==CovariantComponent::Trace         || CComponent==CovariantComponent::Traceless || CComponent==CovariantComponent::Hodge ) trace[n] = _g * ( _gAdj * cov[n] ).trace() / 2.;
				}

				SquareMatrix< double , NumE > mass;
				// The inner product of two maps M,N: TM -> T^*M is
				//		<M,N> = tr( ( gInv * ( gInv * N ).transpose() * g * ( gInv * M ) );
				//            = tr( ( gInv * N.transpose() * gInv * M );
				for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int m=0 ; m<NumE ; m++ )
					if      constexpr( CComponent==CovariantComponent::Symmetric     ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *   sym[n].transpose() * _gAdj *   sym[m] );
					else if constexpr( CComponent==CovariantComponent::AntiSymmetric ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *  aSym[n].transpose() * _gAdj *  aSym[m] );
					else if constexpr( CComponent==CovariantComponent::Trace         ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj * trace[n].transpose() * _gAdj * trace[m] );
					else if constexpr( CComponent==CovariantComponent::Traceless     ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *   sym[n].transpose() * _gAdj *   sym[m] ) - SquareMatrix< double , K >::Trace( _gAdj * trace[n].transpose() * _gAdj * trace[m] );
					else if constexpr( CComponent==CovariantComponent::Hodge         ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj * trace[n].transpose() * _gAdj * trace[m] ) + SquareMatrix< double , K >::Trace( _gAdj *  aSym[n].transpose() * _gAdj *  aSym[m] );
					else                                                               mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *   cov[n].transpose() * _gAdj *   cov[m] );
				return mass / s;
			};
		};
	};
	return this->template system< QuadratureSamples >( SimplexProcessing::ArrayWrapper( Sys ) , false );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , CovariantComponent CComponent , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::stiffness( ScaleFactorField && SF ) const
{
	// Represent the covariant derivative as a map from TM -> T^*M
	static const unsigned int NumE = SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements;

	auto Sys = [&]( size_t sIdx )
		{
			return [tI2E=SimplexProcessing::PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >(simplexVertices(sIdx),simplexNormals(sIdx)),tG=this->metricTensorField(sIdx),tS=this->measureScaleField(sIdx),tSF=SF[sIdx]]( SimplexProcessing::Position< K > p )
				{
					return [i2e=tI2E(p),g=tG(p),s=tS(p),sf=tSF(p)]( const SimplexProcessing::Differential< K , Vector > d[] )
						{
							SquareMatrix< double , K > _gAdj = g.adjugate() / s , _g = g / s;
							SquareMatrix< double , K > cov[ NumE ] , sym[ NumE ] , aSym[ NumE ] , trace[ NumE ];
							Point< double , K > e;
							for( unsigned int n=0 ; n<NumE ; n++ )
							{
								for( unsigned int k=0 ; k<K ; k++ )
								{
									e[k] = 1.;
									Point< double , Dim > _d = i2e * e;
									for( unsigned int _k=0 ; _k<K ; _k++ ) cov[n](k,_k) = Point< double , Dim >::Dot( _d , d[n][_k] );
									e[k] = 0.;
								}

								if constexpr( CComponent==CovariantComponent::AntiSymmetric || CComponent==CovariantComponent::Hodge     )  aSym[n] = ( cov[n] - cov[n].transpose() ) / 2.;
								if constexpr( CComponent==CovariantComponent::Symmetric     || CComponent==CovariantComponent::Traceless )   sym[n] = ( cov[n] + cov[n].transpose() ) / 2.;
								if constexpr( CComponent==CovariantComponent::Trace         || CComponent==CovariantComponent::Traceless || CComponent==CovariantComponent::Hodge ) trace[n] = _g * ( _gAdj * cov[n] ).trace() / 2.;
							}

							SquareMatrix< double , NumE > mass;
							// The inner product of two maps M,N: TM -> T^*M is
							//		<M,N> = tr( ( gInv * ( gInv * N ).transpose() * g * ( gInv * M ) );
							//            = tr( ( gInv * N.transpose() * gInv * M );
							for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int m=0 ; m<NumE ; m++ )
								if      constexpr( CComponent==CovariantComponent::Symmetric     ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *   sym[n].transpose() * _gAdj *   sym[m] );
								else if constexpr( CComponent==CovariantComponent::AntiSymmetric ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *  aSym[n].transpose() * _gAdj *  aSym[m] );
								else if constexpr( CComponent==CovariantComponent::Trace         ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj * trace[n].transpose() * _gAdj * trace[m] );
								else if constexpr( CComponent==CovariantComponent::Traceless     ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *   sym[n].transpose() * _gAdj *   sym[m] ) - SquareMatrix< double , K >::Trace( _gAdj * trace[n].transpose() * _gAdj * trace[m] );
								else if constexpr( CComponent==CovariantComponent::Hodge         ) mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj * trace[n].transpose() * _gAdj * trace[m] ) + SquareMatrix< double , K >::Trace( _gAdj *  aSym[n].transpose() * _gAdj *  aSym[m] );
								else                                                               mass(n,m) = SquareMatrix< double , K >::Trace( _gAdj *   cov[n].transpose() * _gAdj *   cov[m] );
							return mass * sf / s;
						};
				};
		};
	return this->template system< QuadratureSamples >( SimplexProcessing::ArrayWrapper( Sys ) , false );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , typename TangentVectorField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::bracketEnergy( TangentVectorField && V ) const
{
	static const unsigned int NumE = SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements;

	auto Sys = [&]( size_t sIdx )
	{
		return [tE2I=SimplexProcessing::PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >(simplexVertices(sIdx),simplexNormals(sIdx)),tN=this->normalField(sIdx),tG=this->metricTensorField(sIdx),tGInv=this->inverseMetricTensorField(sIdx),tV=V[sIdx]]( SimplexProcessing::Position< K > p )
		{
			return [e2i=tE2I(p),n=tN(p),gInv=tGInv(p),v=tV(p),dv=tV.d(p)]( const std::pair< Vector , SimplexProcessing::Differential< K , Vector > > d[] )
			{
				Point< double , K > _x = e2i * v;
				Point< double , K+1 > _z[ NumE ];
				for( unsigned int e=0 ; e<NumE ; e++ )
				{
					_z[e] = d[e].second( _x ) - dv( e2i * d[e].first );
					_z[e] -= n * Point< double , K+1 >::Dot( _z[e] , n );
				}

				SquareMatrix< double , NumE > mass;
				for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int m=0 ; m<NumE ; m++ )
					mass(n,m) = Point< double , K+1 >::Dot( _z[n] , _z[m] );
				return mass;
			};
		};
	};
	return this->template system< QuadratureSamples >( SimplexProcessing::ArrayWrapper( Sys ) );
}

template< unsigned int K >
template< unsigned int QuadratureSamples , typename TangentVectorField >
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::dotProductEnergy( TangentVectorField && V ) const
{
	static const unsigned int NumE = SimplexProcessing::PhongRodriguesVectorElements< K , Dim >::NumElements;

	auto Sys = [&]( size_t sIdx )
	{
		return [tV=V[sIdx]]( SimplexProcessing::Position< K > p )
		{
			return [v=tV(p)]( const Vector d[] )
			{
				double dot[ NumE ];
				for( unsigned int n=0 ; n<NumE ; n++ ) dot[n] = Point< double , K+1 >::Dot( v , d[n] );

				SquareMatrix< double , NumE > mass;
				for( unsigned int n=0 ; n<NumE ; n++ ) for( unsigned int m=0 ; m<NumE ; m++ ) mass(n,m) = dot[n] * dot[m];
				return mass;
			};
		};
	};
	return this->template system< QuadratureSamples >( SimplexProcessing::ArrayWrapper( Sys ) );
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
Eigen::SparseMatrix< double > EmbeddedPhongMesh< K >::rotation( double radians ) const requires( K==2 )
{
	auto NormalRotation = []( Point< double , Dim > n , double radians )
	{
		double r = cos( radians / 2. );
		Point< double , 3 > i = n * sin( radians / 2. );
		SquareMatrix< double , Dim > R;

		R(0,0) = 1 - 2 * ( i[1]*i[1] + i[2]*i[2] );
		R(1,1) = 1 - 2 * ( i[0]*i[0] + i[2]*i[2] );
		R(2,2) = 1 - 2 * ( i[0]*i[0] + i[1]*i[1] );

		R(1,0) = 2 * ( i[0]*i[1] - r * i[2] );
		R(2,0) = 2 * ( i[0]*i[2] + r * i[1] );
		R(2,1) = 2 * ( i[1]*i[2] - r * i[0] );

		R(0,1) = 2 * ( i[0]*i[1] + r * i[2] );
		R(0,2) = 2 * ( i[0]*i[2] - r * i[1] );
		R(1,2) = 2 * ( i[1]*i[2] + r * i[0] );

		return R;
	};

	Eigen::SparseMatrix< double > Rot;
	std::vector< Eigen::Triplet< double > > triplets( _normals.size()*Dim*Dim );

	ThreadPool::ParallelFor
	(
		0 , _normals.size() ,
		[&]( size_t i )
		{
			SquareMatrix< double , Dim > R = NormalRotation( _normals[i] / Point< double , Dim >::Length( _normals[i] ) , radians );
			for( unsigned int j=0 ; j<Dim ; j++ ) for( unsigned int k=0 ; k<Dim ; k++ )
				triplets[ i*Dim*Dim + j*Dim + k ] = Eigen::Triplet< double >( static_cast< int >( i*Dim+j ) , static_cast< int >( i*Dim+k ) , R(j,k) );
		}
	);

	Rot.resize( _normals.size()*Dim , _normals.size()*Dim );
	Rot.setFromTriplets( triplets.begin() , triplets.end() );
	return Rot;
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
