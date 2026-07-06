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

#ifdef USE_ADJACENCY
///////////////
// Adjacency //
///////////////

size_t Adjacency::oppositeCorner( size_t c ) const
{
	size_t t ; unsigned int k;
	FactorCornerIndex( c , t , k );
	size_t v = cornerToVertex( c );
	size_t e = cornerToEdge( c );
	if     ( cornerToVertex( _dualEdgeCorners[e].first )==v ) return _dualEdgeCorners[e].second;
	else if( cornerToVertex( _dualEdgeCorners[e].second )==v ) return _dualEdgeCorners[e].first;
	else
	{
		MK_THROW( "Not an edge corner" );
		return -1;
	}
}

size_t Adjacency::cornerToVertex( size_t c ) const
{
	size_t t ; unsigned int k;
	FactorCornerIndex( c , t , k );
	return _triangles[t][k];
}

std::vector< size_t > Adjacency::oneRingCorners( size_t v ) const
{
	std::vector< size_t > oneRing;

	size_t c = _vertexToCorner[v];
	size_t t ; unsigned int k;
	FactorCornerIndex( c , t , k );

	size_t nbr = CornerIndex( t , (k+1)%(K+1) );

	oneRing.push_back( nbr );
	while( true )
	{
		size_t nextNbr;

		size_t opposite = oppositeCorner( nbr );
		size_t t ; unsigned int k;
		FactorCornerIndex( opposite , t , k );
		if     ( _triangles[t][(k+1)%(K+1)]==v ) nextNbr = CornerIndex( t , (k+2)%(K+1) );
		else if( _triangles[t][(k+2)%(K+1)]==v ) nextNbr = CornerIndex( t , (k+1)%(K+1) );
		else MK_THROW( "Could not find vertex in neighboring triangle" );
		if( nextNbr==oneRing[0] ) break;
		else oneRing.push_back( nextNbr );

		nbr = nextNbr;
	}

	return oneRing;
}

std::vector< size_t > Adjacency::oneRingVertices( size_t v ) const
{
	std::vector< size_t > oneRing = oneRingCorners( v );
	for( unsigned int i=0 ; i<oneRing.size() ; i++ ) oneRing[i] = cornerToVertex( oneRing[i] );
	return oneRing;
}

Adjacency::Adjacency( const std::vector< SimplexIndex< K > > & triangles , size_t vNum ) : _vNum(vNum) , _triangles(triangles)
{
	auto GetEdge =[]( size_t v1 , size_t v2 ){ return v1<v2 ? std::make_pair( v1 , v2 ) : std::make_pair( v2 , v1 ); };
	struct EdgeData
	{
		EdgeData( size_t idx=-1 , size_t c1=-1 , size_t c2=-1 ) : idx(idx) , c1(c1) , c2(c2){}
		size_t idx , c1 , c2;
	};

	std::map< Edge , EdgeData > edgeMap;

	for( size_t i=0 ; i<triangles.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
	{
		Edge e = GetEdge( _triangles[i][ (k+1)%(K+1) ] , _triangles[i][ (k+2)%(K+1) ] );
		auto it = edgeMap.find(e);
		if( it==edgeMap.end() ) edgeMap[e] = EdgeData( static_cast< unsigned int >( edgeMap.size() ) , CornerIndex( i , k ) );
		else
		{
			if( it->second.c2!=-1 ) MK_THROW( "More than two triangles incident on edge" );
			it->second.c2 = CornerIndex( i , k );
		}
	}
	for( auto it=edgeMap.begin() ; it!=edgeMap.end() ; it++ ) if( it->second.c2==-1 ) MK_THROW( "Boundary edges not supported" );

	_edges.resize( edgeMap.size() );
	_dualEdgeCorners.resize( edgeMap.size() );
	_cornerToEdge.resize( _triangles.size()*(K+1) );

	for( auto it=edgeMap.begin() ; it!=edgeMap.end() ; it++ )
	{
		_edges[ it->second.idx ] = it->first;
		_dualEdgeCorners[ it->second.idx ] = std::make_pair( it->second.c1 , it->second.c2 );
		_cornerToEdge[ it->second.c1 ] = _cornerToEdge[ it->second.c2 ] = it->second.idx;
	}

	_vertexToCorner.resize( _vNum , -1 );
	for( size_t t=0 ; t<_triangles.size() ; t++ ) for( unsigned int k=0 ; k<=K ; k++ ) _vertexToCorner[ _triangles[t][k] ] = CornerIndex( t , k );
}
#endif // USE_ADJACENCY

////////////////
// Subdivider //
////////////////

Eigen::MatrixXd Subdivider::LoopStencil( unsigned int valence , bool useWarren )
{
	auto BlendBeta = [&]( unsigned int v )
	{
		if( useWarren )
		{
			if( v==3 ) return 3./16;
			else       return 3./(8*v);
		}
		else
		{
			double temp = ( 3./8 + 1./4 * cos( 2.* M_PI / v )  );
			return 1./v * ( 5./8 - temp * temp);
		}
	};

	double beta = BlendBeta( valence );

	Eigen::MatrixXd stencil( valence+1 , valence+1 );
	stencil.setZero();
	stencil(0,0) = ( 1. - valence * beta );
	for( unsigned int c=0 ; c<valence ; c++ ) stencil(0,c+1) = beta;
	for( unsigned int r=0 ; r<valence ; r++ )
	{
		stencil(r+1,0) = stencil(r+1,r+1) = 3./8;
		stencil( r+1 , 1+( (r+1)%valence ) ) = stencil( r+1 , 1 + ( (r+valence-1)%valence ) ) = 1./8;
	}
	return stencil;
}


#ifdef USE_ADJACENCY
Subdivider::Subdivider( const std::vector< SimplexIndex< K > > & triangles , size_t vNum ) : _adjacency( triangles , vNum )
{
	_valence.resize( vNum , 0 );
	for( unsigned int i=0 ; i<triangles.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ ) _valence[ triangles[i][k] ]++;
}
#else // !USE_ADJACENCY
Subdivider::Subdivider( const std::vector< SimplexIndex< K > > & triangles , size_t vNum ) : _vNum(vNum) , _triangles(triangles)
{
	auto CornerIndex = []( unsigned int t , unsigned int k ){ return t*(K+1)+k; };
	auto FactorCornerIndex = []( unsigned int c , unsigned int & t , unsigned int & k ){ t = c/(K+1) , k = c%(K+1); };
	auto GetEdge = []( unsigned int v1 , unsigned int v2 ){ return v1<v2 ? std::make_pair( v1 , v2 ) : std::make_pair( v2 , v1 ); };
	auto OppositeCorner = [&]( unsigned int c )
	{
		unsigned int t , k;
		FactorCornerIndex( c , t , k );
		Edge e = GetEdge( triangles[t][(k+1)%(K+1)] , triangles[t][(k+2)%(K+1)] );
		auto iter = _edges.find(e);
		if( iter==_edges.end() ) MK_THROW( "Couldn't find edge" );
		else if( iter->second.c1==c ){ return iter->second.c2; }
		else if( iter->second.c2==c ){ return iter->second.c1; }
		else MK_THROW( "Not an edge corner" );
	};
	auto CornerToVertex = [&]( unsigned int c )
	{
		unsigned int t , k;
		FactorCornerIndex( c , t , k );
		return triangles[t][k];
	};

	_valence.resize( _vNum , 0 );
	for( unsigned int i=0 ; i<triangles.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
	{
		_valence[ triangles[i][k] ]++;

		Edge e = GetEdge( triangles[i][k] , triangles[i][ (k+1)%(K+1) ] );
		auto iter = _edges.find(e);
		if( iter==_edges.end() ) _edges[e] = _EdgeData( static_cast< unsigned int >( _vNum + _edges.size() ) , triangles[i][k] );
		else
		{
			if( iter->second.c2!=-1 ) MK_THROW( "More than two triangles incident on edge" );
			_edges[e].c2 = triangles[i][k];
		}
	}
	for( auto iter=_edges.begin() ; iter!=_edges.end() ; iter++ ) if( iter->second.c2==-1 ) MK_THROW( "Boundary edges not supported" );
}
#endif // USE_ADJACENCY

std::vector< SimplexIndex< 2 > > Subdivider::operator()( void ) const
{
#ifdef USE_ADJACENCY
	std::vector< SimplexIndex< K > > triangles;
	triangles.reserve( _adjacency.tNum() * 4 );
	for( unsigned int i=0 ; i<_adjacency.tNum() ; i++ )
	{
		size_t e[K+1];
		SimplexIndex< K > tri = _adjacency.triangle(i);
		for( unsigned int k=0 ; k<=K ; k++ ) e[k] = _adjacency.vNum() + _adjacency.cornerToEdge( Adjacency::CornerIndex( i , k ) );
		for( unsigned int k=0 ; k<=K ; k++ ) triangles.emplace_back( tri[k] , e[(k+2)%(K+1)] , e[(k+1)%(K+1)] );
		triangles.emplace_back( e[0] , e[1] , e[2] );
	}
	return triangles;
#else // !USE_ADJACENCY
	auto GetEdge = []( unsigned int v1 , unsigned int v2 ){ return v1<v2 ? std::make_pair( v1 , v2 ) : std::make_pair( v2 , v1 ); };

	std::vector< SimplexIndex< K > > simplices;

	simplices.reserve( _triangles.size() * 4 );
	for( unsigned int i=0 ; i<_triangles.size() ; i++ )
	{
		unsigned int e[K+1];
		for( unsigned int k=0 ; k<=K ; k++ )
		{
			auto it = _edges.find( GetEdge( _triangles[i][(k+1)%(K+1)] , _triangles[i][(k+2)%(K+1)] ) );
			if( it==_edges.end() ) MK_THROW( "Edge not in map" );
			e[k] = it->second.idx;
		}
		for( unsigned int k=0 ; k<=K ; k++ ) simplices.emplace_back( _triangles[i][k] , e[(k+2)%(K+1)] , e[(k+1)%(K+1)] );
		simplices.emplace_back( e[0] , e[1] , e[2] );
	}

	return simplices;
#endif // USE_ADJACENCY
}

template< typename T >
std::vector< T > Subdivider::operator()( const std::vector< T > & cData ) const
{
	// Needs to be consistent with the simplex subdivider
	std::vector< T > _cData;
	_cData.reserve( _adjacency.tNum() * 4 * ( K+1 ) );

	for( unsigned int i=0 ; i<_adjacency.tNum() ; i++ )
	{
		T t[K+1] , _t[K+1];
		SimplexIndex< K > tri = _adjacency.triangle(i);
		for( unsigned int k=0 ; k<=K ; k++ ) t[k] = cData[ Adjacency::CornerIndex(i,k) ] , _t[k] = ( cData[ Adjacency::CornerIndex(i,k+1) ] + cData[ Adjacency::CornerIndex(i,k+2) ] ) / 2.;
		for( unsigned int k=0 ; k<=K ; k++ )
		{
			_cData.emplace_back( t[k] );
			_cData.emplace_back( _t[(k+2)%(K+1)] );
			_cData.emplace_back( _t[(k+1)%(K+1)] );
		}
		_cData.emplace_back( _t[0] );
		_cData.emplace_back( _t[1] );
		_cData.emplace_back( _t[2] );
	}
	return _cData;
}

template< typename T >
std::vector< T > Subdivider::operator()( const std::vector< T > & vData , SubdivisionType type ) const
{
#ifdef USE_ADJACENCY
#else // !USE_ADJACENCY
	auto CornerIndex = []( unsigned int t , unsigned int k ){ return t*(K+1)+k; };
	auto FactorCornerIndex = []( unsigned int c , unsigned int & t , unsigned int & k ){ t = c/(K+1) , k = c%(K+1); };
	auto GetEdge = []( unsigned int v1 , unsigned int v2 ){ return v1<v2 ? std::make_pair( v1 , v2 ) : std::make_pair( v2 , v1 ); };
	auto OppositeCorner = [&]( unsigned int c )
	{
		unsigned int t , k;
		FactorCornerIndex( c , t , k );
		Edge e = GetEdge( _triangles[t][(k+1)%(K+1)] , _triangles[t][(k+2)%(K+1)] );
		auto iter = _edges.find(e);
		if( iter==_edges.end() ) MK_THROW( "Couldn't find edge" );
		else if( iter->second.c1==c ){ return iter->second.c2; }
		else if( iter->second.c2==c ){ return iter->second.c1; }
		else MK_THROW( "Not an edge corner" );
	};
	auto CornerToVertex = [&]( unsigned int c )
	{
		unsigned int t , k;
		FactorCornerIndex( c , t , k );
		return _triangles[t][k];
	};
#endif // USE_ADJACENCY

#ifdef USE_ADJACENCY
	if( vData.size()!=_adjacency.vNum() ) MK_THROW( "Unexepcted data size: " , vData.size() , " != " , _adjacency.vNum() );
	std::vector< T > _vData( vData.size() + _adjacency.eNum() );
#else // !USE_ADJACENCY
	if( vData.size()!=_vNum ) MK_THROW( "Unexepcted data size: " , vData.size() , " != " , _vNum );
	std::vector< T > _vData( vData.size() + _edges.size() );
#endif // USE_ADJACENCY
	for( unsigned int i=0 ; i<vData.size() ; i++ ) _vData[i] = vData[i];

#ifdef USE_ADJACENCY
	if( type==SubdivisionType::PLANAR )
		for( size_t e=0 ; e<_adjacency.eNum() ; e++ )
		{
			Adjacency::Edge edge = _adjacency.edge( e );
			_vData[ _adjacency.vNum() + e ] = ( vData[ edge.first ] + vData[ edge.second ] ) / 2;
		}
#else // !USE_ADJACENCY
	if( type==SubdivisionType::PLANAR ) for( auto it=_edges.begin() ; it!=_edges.end() ; it++ ) _vData[ it->second.idx ] = ( vData[ it->first.first ] + vData[ it->first.second ] ) / 2;
#endif // USE_ADJACENCY
	else if( type==SubdivisionType::LOOP_LOOP || type==SubdivisionType::LOOP_WARREN )
	{
#ifdef USE_ADJACENCY
		// Set the edge mid-points
		{
			ThreadPool::ParallelFor
			(
				0 , _adjacency.eNum() ,
				[&]( size_t e )
				{
					Adjacency::Edge edge = _adjacency.edge( e ) , dualEdge = _adjacency.dualEdgeVertices( e );
					_vData[ vData.size()+e ] += 3./8 * ( vData[ edge.first ] + vData[ edge.second ] );
					_vData[ vData.size()+e ] += 1./8 * ( vData[ dualEdge.first ] + vData[ dualEdge.second ] );
				}
			);
		}

		// Update the vertices
		{
			auto BlendBeta = [&]( unsigned int v )
			{
				if( type==SubdivisionType::LOOP_WARREN )
				{
					if( v==3 ) return 3./16;
					else       return 3./(8*v);
				}
				else if( type==SubdivisionType::LOOP_LOOP )
				{
					double temp = ( 3./8 + 1./4 * cos( 2.* M_PI / v )  );
					return 1./v * ( 5./8 - temp * temp);
				}
				else MK_THROW( "Unrecognized subdivision type" );
			};
			ThreadPool::ParallelFor
			(
				0 , _adjacency.vNum() ,
				[&]( size_t v )
				{
#if 0
					std::vector< size_t > oneRing = _adjacency.oneRingCorners( v );
					for( unsigned int i=0 ; i<oneRing.size() ; i++ ) oneRing[i] = _adjacency.cornerToVertex( oneRing[i] );
					Eigen::MatrixXd stencil = LoopStencil( static_cast< unsigned int >( oneRing.size() ) );
					_vData[v] = stencil(0,0) * vData[v];
					for( unsigned int i=0 ; i<oneRing.size() ; i++ ) _vData[v] += stencil(0,1+i) * vData[ oneRing[i] ];
#else
					std::vector< size_t > neighbors = _adjacency.oneRingCorners( v );
					double beta = BlendBeta( static_cast< unsigned int >( neighbors.size() ) );
					_vData[v] = ( 1. - neighbors.size() * beta ) * vData[v];
					for( unsigned int n=0 ; n<neighbors.size() ; n++ ) _vData[v] += beta * vData[ _adjacency.cornerToVertex( neighbors[n] ) ];
#endif
				}
			);
		}


#else // !USE_ADJACENCY
		// Set the edge mid-points
		{
			// Get the contributions from the opposte corner
			for( unsigned int i=0 ; i<_triangles.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
			{
				Edge e = GetEdge( _triangles[i][(k+1)%(K+1)] , _triangles[i][(k+2)%(K+1)] );
				auto it = _edges.find(e);
				if( it==_edges.end() ) MK_THROW( "could not find edge" );
				_vData[ it->second.idx ] += 1./8 * vData[ _triangles[i][k] ];
			}
			// Get the contributions from the end-points
			for( auto it=_edges.begin() ; it!=_edges.end() ; it++ ) _vData[ it->second.idx ] += 3./8 * ( vData[ it->first.first ] + vData[ it->first.second ] );
		}

		// Update the vertices
		{
			auto BlendBeta = [&]( unsigned int v )
			{
				if( type==SubdivisionType::LOOP_WARREN )
				{
					// Using Warren's
					if( v==3 ) return 3./16;
					else       return 3./(8*v);
				}
				else if( type==SubdivisionType::LOOP_LOOP )
				{
					double temp = ( 3./8 + 1./4 * cos( 2.* M_PI / v )  );
					return 1./v * ( 5./8 - temp * temp);
				}
				else MK_THROW( "Unrecognized subdivision type" );
			};

			for( auto iter=_edges.begin() ; iter!=_edges.end() ; iter++ )
			{
				unsigned int v1 = iter->first.first , v2 = iter->first.second;
				double beta1 = BlendBeta( _valence[ v1 ] ) , beta2 = BlendBeta( _valence[ v2 ] );
				_vData[v1] += beta1 * ( vData[v2] - vData[v1] );
				_vData[v2] += beta2 * ( vData[v1] - vData[v2] );
			}
		}
#endif // USE_ADJACENCY
	}
	else if( type==SubdivisionType::BUTTERFLY )
	{
		static const double w = 1./16;
#ifdef USE_ADJACENCY
		for( size_t e=0 ; e<_adjacency.eNum() ; e++ )
		{
			T & v = _vData[ _adjacency.vNum() + 1 ];

			// Edge vertices
			Adjacency::Edge edge = _adjacency.edge( e );
			
			// Opposite vertices
			{
				Adjacency::Edge dualEdge = _adjacency.dualEdgeVertices( e );
				v = ( vData[ edge.first ] + vData[ edge.second ] ) / 2;
				v += ( vData[ dualEdge.first ] + vData[ dualEdge.second ] ) * 2. * w;
			}

			// Flap vertices
			{
				size_t t ; unsigned int k;
				Adjacency::Edge _e = _adjacency.dualEdgeCorners( e );
				size_t c[4];
				{
					Adjacency::FactorCornerIndex( _e.first , t , k );
					c[0] = Adjacency::CornerIndex( t , (k+1)%(K+1) );
					c[1] = Adjacency::CornerIndex( t , (k+1)%(K+2) );
				}
				{
					Adjacency::FactorCornerIndex( _e.second , t , k );
					c[2] = Adjacency::CornerIndex( t , (k+1)%(K+1) );
					c[3] = Adjacency::CornerIndex( t , (k+1)%(K+2) );
				}
				for( unsigned int i=0 ; i<4 ; i++ ) v -= vData[ _adjacency.cornerToVertex( _adjacency.oppositeCorner( c[i] ) ) ] * w;
			}
		}
		MK_ERROR_OUT( "Method not supported" );
#else // !USE_ADJACENCY
		for( auto iter=_edges.begin() ; iter!=_edges.end() ; iter++ )
		{
			// The opposite corners
			unsigned int c[] = { iter->second.c1 , iter->second.c2 };

			// The contribution from the vertices on the edge
			_vData[ iter->second.idx ] = ( vData[ iter->first.first ] + vData[ iter->first.second ] ) / 2.;

			// The contribtion from the opposite corners
			_vData[ iter->second.idx ] += 2. * w * ( vData[ CornerToVertex( c[0] ) ] + vData[ CornerToVertex( c[1] ) ] );

			// The contributions from the flaps
			for( unsigned int i=0 ; i<2 ; i++ )
			{
				// The other two corners in c[i]'s triangle
				unsigned int _c[2];
				{
					unsigned int t , k;
					FactorCornerIndex( c[i] , t , k );
					_c[0] = CornerIndex( t , (k+1+i)%K+1 );
					_c[0] = CornerIndex( t , (k+2+i)%K+1 );
				}

				for( unsigned int j=0 ; j<2 ; j++ )
				{
					unsigned int idx = CornerToVertex( OppositeCorner( _c[j] ) );
					_vData[ iter->second.idx ] -= w * vData[idx];
				}
			}
		}
#endif // USE_ADJACENCY
	}
	return _vData;
}

template< unsigned int Dim >
void Subdivide( std::vector< Point< double , Dim > > & vertices , std::vector< SimplexIndex< 2 > > & simplices , SubdivisionType type )
{
	Subdivider sub( simplices , vertices.size() );

	vertices = sub( vertices , type );

	// [WARNING] This has to be called last because Subdivider stores a reference to the simplices
	simplices = sub();
}

template< unsigned int Dim >
void Subdivide( std::vector< Point< double , Dim > > & vertices , std::vector< Point< double , Dim > > & normals , std::vector< SimplexIndex< 2 > > & simplices , SubdivisionType type )
{
	Subdivider sub( simplices , vertices.size() );

	vertices = sub( vertices , type );
	normals = sub( normals , type );

	// [WARNING] This has to be called last because Subdivider stores a reference to the simplices
	simplices = sub();
}
