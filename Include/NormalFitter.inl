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

template< unsigned int K >
Eigen::MatrixXd NormalFitter< K >::_NormalLimitStencil( Eigen::MatrixXd stencil , double eps )
{
	unsigned int valence = static_cast< unsigned int >( stencil.rows() )-1;
	Eigen::MatrixXd _stencil( valence , valence );
	for( unsigned int i=0 ; i<valence ; i++ ) for( unsigned int j=0 ; j<valence ; j++ ) _stencil(i,j) = stencil(i+1,j+1) - stencil(0,j+1);
	Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eSolver( _stencil );
	Eigen::VectorXd eValues = eSolver.eigenvalues();
	for( unsigned int i=0 ; i<eValues.size() ; i++ ) eValues[i] = ( eValues[valence-1]-eValues[i] )<eps ? 1. : 0.;
	return eSolver.eigenvectors() * eValues.asDiagonal() * eSolver.eigenvectors().transpose();
}

template< unsigned int K >
std::vector< Point< double , K+1 > > NormalFitter< K >::_NormalsFromSimplices( const std::vector< Point< double , Dim > > & vertices , const std::vector< SimplexIndex< K > > & simplices , unsigned int iters , bool useWarren )
{
	std::map< unsigned int , Eigen::MatrixXd > valenceToStencil;
	Subdivide::Adjacency adjacency( simplices , vertices.size() );
	std::vector< Point< double , Dim > > normals( vertices.size() );

	for( unsigned int v=0 ; v<vertices.size() ; v++ )
	{
		std::vector< size_t > oneRing = adjacency.oneRingVertices( v );
		unsigned int valence = static_cast< unsigned int >( oneRing.size() );
		Eigen::MatrixXd stencil = Eigen::MatrixXd(valence+1,valence+1);
		stencil.setIdentity();

		if( iters )
		{
			auto it = valenceToStencil.find( valence );
			if( it==valenceToStencil.end() )
			{
				Eigen::MatrixXd _stencil = Subdivide::Subdivider::LoopStencil( valence , useWarren )*2;
				for( unsigned int i=0 ; i<iters ; i++ ) stencil = stencil * _stencil;
				valenceToStencil[ valence ] = stencil;
			}
			else stencil = it->second;
		}

		std::vector< Point< double , 3 > > limitVertices( oneRing.size()+1 );

		limitVertices[0] = vertices[v] * stencil(0,0);
		for( unsigned int c=0 ; c<oneRing.size() ; c++ ) limitVertices[0] += vertices[ oneRing[c] ] * stencil(0,c+1);
		for( unsigned int r=0 ; r<oneRing.size() ; r++ )
		{
			limitVertices[r+1] = vertices[v] * stencil(r+1,0);
			for( unsigned int c=0 ; c<oneRing.size() ; c++ ) limitVertices[r+1] += vertices[ oneRing[c] ] * stencil(r+1,c+1);
		}

		Simplex< double , Dim , K > s;
		Point< double , Dim > n;
		s[0] = limitVertices[0];
		for( unsigned int i=0 ; i<oneRing.size() ; i++ )
		{
			s[1] = limitVertices[1+i];
			s[2] = limitVertices[1+(i+1)%oneRing.size() ];
			n += s.normal();
		}
		normals[v] = n;
	}
	return normals;
}

template< unsigned int K >
std::vector< Point< double , K+1 > > NormalFitter< K >::_LimitNormalsFromSimplices( const std::vector< Point< double , Dim > > & vertices , const std::vector< SimplexIndex< K > > & simplices , bool useWarren , double eps )
{
	std::map< unsigned int , Eigen::MatrixXd > valenceToStencil;
	Subdivide::Adjacency adjacency( simplices , vertices.size() );
	std::vector< Point< double , Dim > > normals( vertices.size() );

	for( unsigned int v=0 ; v<vertices.size() ; v++ )
	{
		std::vector< size_t > oneRing = adjacency.oneRingVertices( v );
		unsigned int valence = static_cast< unsigned int >( oneRing.size() );
		Eigen::MatrixXd stencil = Eigen::MatrixXd(valence,valence);
		stencil.setIdentity();

		auto it = valenceToStencil.find( valence );
		if( it==valenceToStencil.end() ) valenceToStencil[ valence ] = stencil = _NormalLimitStencil( Subdivide::Subdivider::LoopStencil( valence , useWarren ) , eps );
		else stencil = it->second;

		std::vector< Point< double , 3 > > limitVertices( oneRing.size() );
		for( unsigned int r=0 ; r<oneRing.size() ; r++ ) for( unsigned int c=0 ; c<oneRing.size() ; c++ ) limitVertices[r] += ( vertices[ oneRing[c] ] - vertices[v] ) * stencil(r,c);

		Simplex< double , Dim , K > s;
		Point< double , Dim > n;
		for( unsigned int i=0 ; i<oneRing.size() ; i++ )
		{
			s[1] = limitVertices[i];
			s[2] = limitVertices[(i+1)%oneRing.size() ];
			n += s.normal();
		}

		normals[v] = n;
	}
	return normals;
}

template< unsigned int K >
template< unsigned int Quadrature >
void NormalFitter< K >::_SmoothNormals( EmbeddedMesh< K , Dim > mesh , std::vector< Point< double , Dim > > & normals , double diffusionTime )
{
	auto Extend = []( const Eigen::SparseMatrix< double > & M )
		{
			std::vector< Eigen::Triplet< double > > triplets;
			triplets.reserve( M.nonZeros()*Dim );
			for( unsigned int i=0; i<M.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter( M , i ) ; iter ; ++iter )
				for( unsigned int d=0 ; d<Dim ; d++ ) triplets.emplace_back( static_cast< int >( iter.row()*Dim+d ) , static_cast< int >( iter.col()*Dim+d ) , iter.value() );
			Eigen::SparseMatrix< double > L( M.rows() * Dim , M.cols()*Dim );
			L.setFromTriplets( triplets.begin() , triplets.end() );
			return L;
		};

	////////////////
	// [Energy 1] //
	////////////////
	// E(X) = ( P*X )^T * M * ( P*X ) + eps * ( Y + P*X )^T * S * ( Y + P*X )
	//      = X^t * P^t * M * P * X + eps * ( X^t * P^t * S * P * X + 2 * X^t * P^t * S * Y ) + ...
	// Differentiating:
	//    0 = P^t * M * P * X + eps * ( P^t * S * P * X + P^t * S * Y )
	//      = P^t * ( M + eps * S ) * P * X + eps * P^t * S * Y
	//  =>
	//   P^t * ( M + eps * S ) * P * X = - eps * P^t * S * Y

	////////////////
	// [Energy 2] //
	////////////////
	// E(X) = ( Y + P*X )^T * M * ( Y + P*X ) + eps * ( Y + P*X )^T * S * ( Y + P*X )
	//      = X^t * P^t * M * P * X + 2 * X^t * P^t * M * Y + eps * ( X^t * P^t * S * P * X + 2 * X^t * P^t * S * Y ) + ...
	// Differentiating:
	//    0 = P^t * M * P * X + P^t * M * Y + eps * ( P^t * S * P * X + P^t * S * Y )
	//      = P^t * ( M + eps * S ) * P * X + P^t * ( M + eps * S ) * Y
	//  =>
	//   P^t * ( M + eps * S ) * P * X = - P^t * ( M + eps * S ) * Y

	// Going with [Energy 1]

	Eigen::SparseMatrix< double > L , P , Pt , S = mesh.template scalarStiffness< Quadrature >() * diffusionTime;

	// Set the system
	{
		L = Extend( mesh.template scalarMass< Quadrature >() + S );

		EmbeddedPhongMesh< K > _mesh( mesh.vertices() , normals , mesh.simplices() );
		P = _mesh.tangentProlongation();
		Pt = P.transpose();
	}

	// Perform the factorization
	LLtSolver solver( Pt * L * P );
	if( solver.info()!=Eigen::Success ) MK_ERROR_OUT( "Failed to factorize system" );

	Eigen::VectorXd nDual( normals.size()  * Dim );
	ThreadPool::ParallelFor
		(
			0 , normals.size() ,
			[&]( size_t i )
			{
				Point< double , Dim > dual;
				for( Eigen::SparseMatrix< double >::InnerIterator iter( S , i ) ; iter ; ++iter ) dual += iter.value() * normals[ iter.row() ];
				for( unsigned int d=0 ; d<Dim ; d++ ) nDual[i*Dim+d] = dual[d];
			}
		);

	Eigen::VectorXd v = P * solver.solve( - Pt * nDual );
	for( unsigned int i=0 ; i<normals.size() ; i++ )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) normals[i][d] += v[i*Dim+d];
		normals[i] /= Point< double , Dim >::Length( normals[i] );
	}
}

template< unsigned int K >
template< unsigned int Quadrature , unsigned int ... Quadratures >
void NormalFitter< K >::_SmoothNormals( unsigned int quadrature , EmbeddedMesh< K , Dim > mesh , std::vector< Point< double , Dim > > & normals , double diffusionTime )
{
	if( quadrature==Quadrature ) return _SmoothNormals< Quadrature >( mesh , normals , diffusionTime );
	else if constexpr( sizeof...(Quadratures) ) return _SmoothNormals< Quadratures... >( quadrature , mesh , normals , diffusionTime );
	else MK_THROW( "Bad quadrature option: " , quadrature );
}

template< unsigned int K >
std::vector< Point< double , K+1 > > NormalFitter< K >::Fit( std::vector< Point< double , Dim > > & vertices , std::vector< SimplexIndex< K > > & simplices , Params params )
{
	std::vector< Point< double , Dim > > normals;

	if( params.subdivisionIterations==-1 ) normals = _LimitNormalsFromSimplices( vertices , simplices , params.useWarren , params.spectralEps );
	else                                   normals = _NormalsFromSimplices( vertices , simplices , params.subdivisionIterations , params.useWarren );
	for( unsigned int i=0 ; i<normals.size() ; i++ ) normals[i] /= Point< double , Dim >::Length( normals[i] );

	if( params.diffusionTime )
	{
		if constexpr( K==2 )
		{
			EmbeddedMesh< K , Dim > mesh( vertices , simplices );

			auto AverageEdgeLength = [&]( void )
				{
					size_t count = 0;
					double avg = 0;
					for( unsigned int i=0 ; i<simplices.size() ; i++ ) for( unsigned int j=0 ; j<=K ; j++ ) for( unsigned int k=0 ; k<j ; k++ )
					{
						avg += Point< double , Dim >::Length( vertices[ simplices[i][j] ] - vertices[ simplices[i][k] ] );
						count++;
					}
					return avg / count;
				};
			if( params.normalizeScale )
			{
				MeshNormalizer< Dim > normalizer( mesh , vertices );
				_SmoothNormals< SUPPORTED_QUADRATURE >( params.quadrature , mesh , normals , params.scaleByAverageEdgeLength ? params.diffusionTime * AverageEdgeLength() : params.diffusionTime );
			}
			else _SmoothNormals< SUPPORTED_QUADRATURE >( params.quadrature , mesh , normals , params.scaleByAverageEdgeLength ? params.diffusionTime * AverageEdgeLength() : params.diffusionTime );
		}
		else MK_WARN( "Quadrature only supported for dimension K=2" );
	}
	return normals;
}

