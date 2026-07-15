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

////////////////
// _MeshField //
////////////////

template< unsigned int K , bool NeedsVertices , bool NeedsNormals , typename SimplexField , typename ... InTypes >
_MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... >::_MeshField( const EmbeddedMesh< K , Dim > & mesh , const std::vector< InTypes > & ... in ) requires ( !NeedsNormals )
	: _simplices( mesh.simplices() ) , _vertices( mesh.vertices() ) , _in( std::tie( in ... ) )
{}

template< unsigned int K , bool NeedsVertices , bool NeedsNormals , typename SimplexField , typename ... InTypes >
_MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... >::_MeshField( const EmbeddedPhongMesh< K > & mesh , const std::vector< InTypes > & ... in )
	: _simplices( mesh.simplices() ) , _vertices( mesh.vertices() ) , _normals( mesh.normals() ) , _in( std::tie( in ... ) )
{}

template< unsigned int K , bool NeedsVertices , bool NeedsNormals , typename SimplexField , typename ... InTypes >
auto _MeshField< K , NeedsVertices , NeedsNormals , SimplexField , InTypes ... >::operator[]( size_t sIdx ) const
{
	std::tuple< std::array< InTypes , K+1 > ... > in;
	if      constexpr( sizeof...(InTypes)==1 ) for( unsigned int k=0 ; k<=K ; k++ )
	{
		std::get< 0 >( in )[k] = std::get< 0 >( _in )[ _simplices[sIdx][k] ];
	}
	else if constexpr( sizeof...(InTypes)==2 ) for( unsigned int k=0 ; k<=K ; k++ )
	{
		std::get< 0 >( in )[k] = std::get< 0 >( _in )[ _simplices[sIdx][k] ];
		std::get< 1 >( in )[k] = std::get< 1 >( _in )[ _simplices[sIdx][k] ];
	}

	if constexpr( NeedsVertices )
	{
		if constexpr( NeedsNormals )
		{
			Simplex< double , Dim , K > v , n;
			for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ] , n[k] = _normals[ _simplices[sIdx][k] ];
			if      constexpr( sizeof...(InTypes)==0 ) return SimplexField( v , n );
			else if constexpr( sizeof...(InTypes)==1 ) return SimplexField( v , n , &( std::get< 0 >( in ) )[0] );
			else if constexpr( sizeof...(InTypes)==2 ) return SimplexField( v , n , &( std::get< 0 >( in ) )[0] , &( std::get< 1 >( in ) )[0] );
		}
		else
		{
			Simplex< double , Dim , K > v;
			for( unsigned int k=0 ; k<=K ; k++ ) v[k] = _vertices[ _simplices[sIdx][k] ];
			if      constexpr( sizeof...(InTypes)==0 ) return SimplexField( v );
			else if constexpr( sizeof...(InTypes)==1 ) return SimplexField( v , &( std::get< 0 >( in ) )[0] );
			else if constexpr( sizeof...(InTypes)==2 ) return SimplexField( v , &( std::get< 0 >( in ) )[0] , &( std::get< 1 >( in ) )[0] );
		}
	}
	else
	{
		if constexpr( NeedsNormals )
		{
			Simplex< double , Dim , K > n;
			for( unsigned int k=0 ; k<=K ; k++ ) n[k] = _normals[ _simplices[sIdx][k] ];
			if      constexpr( sizeof...(InTypes)==0 ) return SimplexField( n );
			else if constexpr( sizeof...(InTypes)==1 ) return SimplexField( n , &( std::get< 0 >( in ) )[0] );
			else if constexpr( sizeof...(InTypes)==2 ) return SimplexField( n , &( std::get< 0 >( in ) )[0] , &( std::get< 1 >( in ) )[0] );
		}
		else
		{
			if      constexpr( sizeof...(InTypes)==0 ) return SimplexField();
			else if constexpr( sizeof...(InTypes)==1 ) return SimplexField( &( std::get< 0 >( in ) )[0] );
			else if constexpr( sizeof...(InTypes)==2 ) return SimplexField( &( std::get< 0 >( in ) )[0] , &( std::get< 1 >( in ) )[0] );
		}
	}
}

template< unsigned int K , SimplexProcessing::HasDotProduct T , HasMeshFunction< K , T > Field >
auto SquaredL2NormField( Field && F )
{
	return SimplexProcessing::ArrayWrapper( [ f=std::forward< Field >(F) ]( size_t sIdx ){ return SimplexProcessing::SquaredL2NormField< K , T >( f[sIdx] ); } );
}

template< unsigned int K , SimplexProcessing::HasDotProduct T , HasMeshFunction< K , T > Field1 , HasMeshFunction< K , T > Field2 >
auto SquaredL2DifferenceField( Field1 && F1 , Field2 && F2 )
{
	return SimplexProcessing::ArrayWrapper( [ f1=std::forward< Field1 >(F1) , f2=std::forward< Field2 >(F2) ]( size_t sIdx ){ return SimplexProcessing::SquaredL2DifferenceField< K , T  >( f1[sIdx] , f2[sIdx] ); } );
}

template< unsigned int K , SimplexProcessing::HasDotProduct T , HasMeshFunction< K , T > Field1 , HasMeshFunction< K , T > Field2 >
auto DotProductField( Field1 && F1 , Field2 && F2 )
{
	return SimplexProcessing::ArrayWrapper( [ f1=std::forward< Field1 >(F1) , f2=std::forward< Field2 >(F2) ]( size_t sIdx ){ return SimplexProcessing::DotProductField< K , T >( f1[sIdx] , f2[sIdx] ); } );
}

template< unsigned int K , HasMeshFunction< K , Point< double , K > > VectorField , HasMeshDifferentiableFunction< K , double > ScalarField >
auto DerivationField( ScalarField && SF , VectorField && VF )
{
	return SimplexProcessing::ArrayWrapper( [ sf=std::forward< ScalarField >(SF) , vf=std::forward< VectorField >(VF) ]( size_t sIdx ){ return SimplexProcessing::DerivationField< K >( sf[sIdx] , vf[sIdx] ); } );
}