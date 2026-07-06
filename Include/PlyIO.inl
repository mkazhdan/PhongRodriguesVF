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

////////////////////////
// Read Functionality //
////////////////////////
template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void _ReadMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< SimplexIndex< K , Index > > &simplices , int & file_type )
{
	using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > >;
	using Vertex = typename VertexFactory::DataType;
	VertexFactory factory;

	std::vector< Vertex > _vertices;
	std::vector< bool > readFlags;
	file_type = PLY::ReadSimplices( fileName , factory , _vertices , simplices , readFlags );

	vertices.resize( _vertices.size() );
	for( unsigned int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = _vertices[i].template get<0>();
	if( factory.template plyValidReadProperties< 1 >( readFlags ) )
	{
		normals.resize( _vertices.size() );
		for( unsigned int i=0 ; i<_vertices.size() ; i++ ) normals[i] = _vertices[i].template get<1>();
	}
	else if constexpr( Dim==K+1 ) normals = SimplexToVertexNormals( vertices , simplices );
	else MK_THROW( "Could not compute normals for non-co-dimension-one manifold: " , Dim , " != " , K , "+1" );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void _ReadColorMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< Point< Real , 3 > > &colors , std::vector< SimplexIndex< K , Index > > &simplices , int & file_type )
{
	using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > >;
	using Vertex = typename VertexFactory::DataType;
	VertexFactory factory;

	std::vector< Vertex > _vertices;
	std::vector< bool > readFlags;
	file_type = PLY::ReadSimplices( fileName , factory , _vertices , simplices , readFlags );

	vertices.resize( _vertices.size() ) , colors.resize( _vertices.size() );
	for( unsigned int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = _vertices[i].template get<0>() , colors[i] = _vertices[i].template get<2>() / static_cast< Real >( 255. );
	if( factory.template plyValidReadProperties< 1 >( readFlags ) )
	{
		normals.resize( _vertices.size() );
		for( unsigned int i=0 ; i<_vertices.size() ; i++ ) normals[i] = _vertices[i].template get<1>();
	}
	else if constexpr( Dim==K+1 ) normals = SimplexToVertexNormals( vertices , simplices );
	else MK_THROW( "Could not compute normals for non-co-dimension-one manifold: " , Dim , " != " , K , "+1" );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim , unsigned int TDim >
void _ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< K , Index > > &simplices , int & file_type )
{
	bool hasNormals = false;
	if constexpr( K==2 )
	{
		using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > >;
		using Vertex = typename VertexFactory::DataType;
		using Face = PlyTexturedFace< unsigned int , Real >;

		VertexFactory factory;

		std::vector< Vertex > inVertices;
		std::vector< Face > inFaces;
		std::vector< bool > vFlags , fFlags;

		try{ file_type = PLY::ReadPolygons( fileName , factory , inVertices , inFaces , Face::ReadProperties , Face::ReadComponents , vFlags , fFlags ); }
		catch( const MishaK::Exception & ){ MK_ERROR_OUT( "Failed to read ply file: " , fileName ); }

		hasNormals = factory.template plyValidReadProperties< 1 >( vFlags );
		if( !fFlags[0] ) MK_ERROR_OUT( "Failed to read face indices" );
		if( !fFlags[1] ) MK_ERROR_OUT( "Failed to read face textures" );

		size_t faceNum = inFaces.size();
		for( unsigned int i=(unsigned int)inFaces.size() ; i!=0 ; i-- )
		{
			Face &face = inFaces[i-1];
			Face oldFace = face;

			if( face.size()>3 )
			{
				std::vector< Point< Real , Dim > > _vertices( face.size() );
				std::vector< SimplexIndex< K > > _triangles;
				for( unsigned int j=0 ; j<(unsigned int)face.size() ; j++ ) _vertices[j] = inVertices[ face[j] ].template get<0>();
				MinimalAreaTriangulation::GetTriangulation( _vertices , _triangles );

				auto TriangleToFace = [&]( SimplexIndex< K > si )
					{
						Face face;
						face.resize(K+1);
						for( unsigned int k=0 ; k<=K ; k++ ) face[k] = oldFace[ si[k] ] , face.texture(k) = oldFace.texture( si[k] );
						return face;
					};

				face = TriangleToFace( _triangles[0] );
				for( unsigned int j=1 ; j<_triangles.size() ; j++ ) inFaces.push_back( TriangleToFace( _triangles[j] ) );
			}
		}

		for( unsigned int i=0 ; i<inFaces.size() ; i++ )
			if( inFaces[i].nr_uv_coordinates!=(K+1)*K ) MK_ERROR_OUT( "Unexpected number of texture coordinates: " , inFaces[i].nr_uv_coordinates , " != " , K+1 , " * " , K );

		simplices.resize( inFaces.size() );
		textureCoordinates.resize( inFaces.size()*(K+1) );
		for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
		{
			simplices[i][k] = inFaces[i][k];
			textureCoordinates[ i*(K+1)+k ] = inFaces[i].texture(k);
		}

		vertices.resize( inVertices.size() );
		for( unsigned int i=0 ; i<inVertices.size() ; i++ ) vertices[i] = inVertices[i].template get<0>();
		if( hasNormals )
		{
			normals.resize( inVertices.size() );
			for( unsigned int i=0 ; i<inVertices.size() ; i++ ) normals[i] = inVertices[i].template get<1>();
		}
		else if constexpr( K+1==Dim ) normals = SimplexToVertexNormals( vertices , simplices );
	}
	else
	{
	}
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< K , Index > > &simplices )
{
	int file_type;
	try
	{
		std::vector< Point< Real , Dim > > normals;
		_ReadMesh( fileName , vertices , normals , simplices , file_type );
	}
	catch( ... ){}
	return file_type;
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< SimplexIndex< K , Index > > &simplices )
{
	int file_type;
	_ReadMesh( fileName , vertices , normals , simplices , file_type );
	return file_type;
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadColorMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , 3 > > &colors , std::vector< SimplexIndex< K , Index > > &simplices )
{
	int file_type;
	try
	{
		std::vector< Point< Real , Dim > > normals;
		_ReadColorMesh( fileName , vertices , normals , colors , simplices , file_type );
	}
	catch( ... ){}
	return file_type;
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadColorMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< Point< Real , 3 > > &colors , std::vector< SimplexIndex< K , Index > > &simplices )
{
	int file_type;
	_ReadColorMesh( fileName , vertices , normals , colors , simplices , file_type );
	return file_type;
}

template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
int ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< 2 , Index > > &simplices )
{
	static const unsigned int K=2;
	using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > >;
	using Vertex = typename VertexFactory::DataType;
	using Face = PlyTexturedFace< unsigned int , Real >;

	VertexFactory factory;

	std::vector< Vertex > inVertices;
	std::vector< Face > inFaces;
	bool *vFlags = new bool[ factory.plyReadNum() ];
	bool *fFlags = new bool[ Face::ReadComponents ];
	int file_type;

	try{ file_type = PLY::ReadPolygons( fileName , factory , inVertices , inFaces , Face::ReadProperties , Face::ReadComponents , vFlags , fFlags ); }
	catch( const MishaK::Exception & ){ MK_ERROR_OUT( "Failed to read ply file: " , fileName ); }

	bool hasNormals = factory.template plyValidReadProperties< 1 >( vFlags );
	if( !fFlags[0] ) MK_ERROR_OUT( "Failed to read face indices" );
	if( !fFlags[1] ) MK_ERROR_OUT( "Failed to read face textures" );
	delete[] vFlags;
	delete[] fFlags;

	if constexpr( K==2 )
	{
		size_t faceNum = inFaces.size();
		for( unsigned int i=(unsigned int)inFaces.size() ; i!=0 ; i-- )
		{
			Face &face = inFaces[i-1];
			Face oldFace = face;

			if( face.size()>3 )
			{
				std::vector< Point< Real , Dim > > _vertices( face.size() );
				std::vector< SimplexIndex< K > > _triangles;
				for( unsigned int j=0 ; j<(unsigned int)face.size() ; j++ ) _vertices[j] = inVertices[ face[j] ].template get<0>();
				MinimalAreaTriangulation::GetTriangulation( _vertices , _triangles );

				auto TriangleToFace = [&]( SimplexIndex< K > si )
					{
						Face face;
						face.resize(K+1);
						for( unsigned int k=0 ; k<=K ; k++ ) face[k] = oldFace[ si[k] ] , face.texture(k) = oldFace.texture( si[k] );
						return face;
					};

				face = TriangleToFace( _triangles[0] );
				for( unsigned int j=1 ; j<_triangles.size() ; j++ ) inFaces.push_back( TriangleToFace( _triangles[j] ) );
			}
		}
		//		if( inFaces.size()!=faceNum ) MK_WARN( "Triangulated: " , faceNum , " -> " , inFaces.size() );
	}

	for( unsigned int i=0 ; i<inFaces.size() ; i++ )
		if( inFaces[i].nr_vertices!=(K+1) ) MK_ERROR_OUT( "Face is not a simplex" );
		else if( inFaces[i].nr_uv_coordinates!=(K+1)*K ) MK_ERROR_OUT( "Unexpected number of texture coordinates: " , inFaces[i].nr_uv_coordinates , " != " , K+1 , " * " , K );

	simplices.resize( inFaces.size() );
	textureCoordinates.resize( inFaces.size()*(K+1) );
	for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
	{
		simplices[i][k] = inFaces[i][k];
		textureCoordinates[ i*(K+1)+k ] = inFaces[i].texture(k);
	}

	vertices.resize( inVertices.size() );
	for( unsigned int i=0 ; i<inVertices.size() ; i++ ) vertices[i] = inVertices[i].template get<0>();
	if( hasNormals )
	{
		normals.resize( inVertices.size() );
		for( unsigned int i=0 ; i<inVertices.size() ; i++ ) normals[i] = inVertices[i].template get<1>();
	}
	else if constexpr( K+1==Dim ) normals = SimplexToVertexNormals( vertices , simplices );
	return file_type;
}

template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
int ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< 2 , Index > > &simplices )
{
	std::vector< Point< Real , Dim > > normals;
	return ReadTexturedMesh( fileName , vertices , normals , textureCoordinates , simplices );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadVectorFieldMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< Point< Real , Dim > > &vectorField , std::vector< SimplexIndex< K , Index > > &simplices )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > , DataFactory::VectorFieldFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	Factory factory;

	std::vector< Vertex > _vertices;
	std::vector< bool > readFlags;
	int file_type = PLY::ReadSimplices( fileName , factory , _vertices , simplices , readFlags );

	vertices.resize( _vertices.size() );
	for( unsigned int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = _vertices[i].template get<0>();

	if( factory.template plyValidReadProperties< 1 >( readFlags ) )
	{
		normals.resize( _vertices.size() );
		for( unsigned int i=0 ; i<_vertices.size() ; i++ ) normals[i] = _vertices[i].template get<1>();
	}

	if( factory.template plyValidReadProperties< 2 >( readFlags ) )
	{
		vectorField.resize( _vertices.size() );
		for( unsigned int i=0 ; i<_vertices.size() ; i++ ) vectorField[i] = _vertices[i].template get<2>();
	}
	return file_type;
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadVectorFieldMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &vectorField , std::vector< SimplexIndex< K , Index > > &simplices )
{
	std::vector< Point< Real , Dim > > normals;
	return ReadVectorFieldMesh( fileName , vertices , normals , vectorField , simplices );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadVectorFieldMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , Dim > > &normals , std::vector< SimplexIndex< K , Index > > &simplices , std::vector< Point< Real , Dim > > &vectorField )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	Factory factory;

	std::vector< GregTurk::PlyProperty > properties = PlyVFFace< Index , Real , Dim >::ReadProperties();
	std::vector< PlyVFFace< Index , Real , Dim > > polygons;
	std::vector< Vertex > _vertices;
	std::vector< bool > vertexReadFlags , polygonReadFlags;
	int file_type = PLY::ReadPolygons( fileName , factory , _vertices , polygons , &properties[0] , static_cast< int >( properties.size() ) , vertexReadFlags , polygonReadFlags );


	vertices.resize( _vertices.size() );
	for( unsigned int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = _vertices[i].template get<0>();

	if( factory.template plyValidReadProperties< 1 >( vertexReadFlags ) )
	{
		normals.resize( _vertices.size() );
		for( unsigned int i=0 ; i<_vertices.size() ; i++ ) normals[i] = _vertices[i].template get<1>();
	}

	simplices.resize( polygons.size() );
	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		if( polygons[i].size()!=(K+1) ) MK_THROW( "Non-simplicial face: "  , polygons[i].size() , " != " , K+1 );
		for( unsigned int k=0 ; k<=K ; k++ ) simplices[i][k] = polygons[i][k];
	}

	if( ( polygonReadFlags[1] || polygonReadFlags[4] ) && ( polygonReadFlags[2] || polygonReadFlags[5] ) && ( polygonReadFlags[2] || polygonReadFlags[6] ) )
	{
		vectorField.resize( polygons.size() );
		for( unsigned int i=0 ; i<polygons.size() ; i++ ) vectorField[i] = polygons[i].v;
	}

	return file_type;
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
int ReadVectorFieldMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< K , Index > > &simplices , std::vector< Point< Real , Dim > > &vectorField )
{
	std::vector< Point< Real , Dim > > normals;
	return ReadVectorFieldMesh( fileName , vertices , normals , simplices , vectorField );
}


/////////////////////////
// Write Functionality //
/////////////////////////
template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using VertexFactory = DataFactory::PositionFactory< Real , Dim >;
	using Vertex = typename VertexFactory::DataType;
	VertexFactory factory;

	PLY::WriteSimplices( fileName , factory , vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , Dim > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > >;
	using Vertex = typename VertexFactory::DataType;

	if( vertices.size()!=normals.size() ) MK_THROW( "Number of vertices and normals don't match: " , vertices.size() , " != " , normals.size() );

	VertexFactory factory;

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ )
	{
		_vertices[i].template get<0>() = vertices[i];
		_vertices[i].template get<1>() = normals[i];
	}

	PLY::WriteSimplices( fileName , factory , _vertices , simplices , file_type );
}


template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteColorMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , 3 > > &colors , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > >;
	using Vertex = typename VertexFactory::DataType;

	if( vertices.size()!=colors.size() ) MK_THROW( "Number of vertices and colors don't match: " , vertices.size() , " != " , colors.size() );

	VertexFactory factory;

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ )
	{
		_vertices[i].template get<0>() = vertices[i];
		_vertices[i].template get<1>() = colors[i] * 255.;
	}

	PLY::WriteSimplices( fileName , factory , _vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteColorMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , Dim > > &normals , const std::vector< Point< Real , 3 > > &colors , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using VertexFactory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > >;
	using Vertex = typename VertexFactory::DataType;

	if( vertices.size()!=colors.size() ) MK_THROW( "Number of vertices and colors don't match: " , vertices.size() , " != " , colors.size() );
	if( vertices.size()!=normals.size() ) MK_THROW( "Number of vertices and normals don't match: " , vertices.size() , " != " , normals.size() );

	VertexFactory factory;

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ )
	{
		_vertices[i].template get<0>() = vertices[i];
		_vertices[i].template get<1>() = normals[i];
		_vertices[i].template get<2>() = colors[i] * 255.;
	}

	PLY::WriteSimplices( fileName , factory , _vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , Dim > > &vf , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::VectorFieldFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	if( vertices.size()!=vf.size() ) MK_THROW( "Number of vertices and vectors don't match: " , vertices.size() , " != " , vf.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = vf[i];
	PLY::WriteSimplices( fileName , Factory() , _vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , Dim > > &vf , const std::vector< Point< Real , Dim > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > , DataFactory::VectorFieldFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	if( vertices.size()!=vf.size() ) MK_THROW( "Number of vertices and vectors don't match: " , vertices.size() , " != " , vf.size() );
	if( vertices.size()!=normals.size() ) MK_THROW( "Number of vertices and normals don't match: " , vertices.size() , " != " , normals.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = normals[i] , _vertices[i].template get<2>() = vf[i];
	PLY::WriteSimplices( fileName , Factory() , _vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteColorVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , 3 > > &colors , const std::vector< Point< Real , Dim > > &vf , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > , DataFactory::VectorFieldFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	if( vertices.size()!=vf.size() ) MK_THROW( "Number of vertices and vectors don't match: " , vertices.size() , " != " , vf.size() );
	if( vertices.size()!=colors.size() ) MK_THROW( "Number of vertices and colors don't match: " , vertices.size() , " != " , colors.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = colors[i]*255 , _vertices[i].template get<2>() = vf[i];
	PLY::WriteSimplices( fileName , Factory() , _vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteColorVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , 3 > > &colors , const std::vector< Point< Real , Dim > > &vf , const std::vector< Point< Real , Dim > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > , DataFactory::NormalFactory< Real , Dim > , DataFactory::VectorFieldFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	if( vertices.size()!=vf.size() ) MK_THROW( "Number of vertices and vectors don't match: " , vertices.size() , " != " , vf.size() );
	if( vertices.size()!=colors.size() ) MK_THROW( "Number of vertices and colors don't match: " , vertices.size() , " != " , colors.size() );
	if( vertices.size()!=normals.size() ) MK_THROW( "Number of vertices and normals don't match: " , vertices.size() , " != " , normals.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = colors[i]*255 , _vertices[i].template get<2>() = normals[i] , _vertices[i].template get<3>() = vf[i];
	PLY::WriteSimplices( fileName , Factory() , _vertices , simplices , file_type );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< SimplexIndex< K , Index > > &simplices , const std::vector< Point< Real , Dim > > &vf , int file_type )
{
	using Factory = DataFactory::PositionFactory< Real , Dim >;
	using Vertex = typename Factory::DataType;

	if( simplices.size()!=vf.size() ) MK_THROW( "Number of simplices and vectors don't match: " , simplices.size() , " != " , vf.size() );

	std::vector< GregTurk::PlyProperty > properties = PlyVFFace< Index , Real , Dim >::Properties();
	std::vector< PlyVFFace< Index , Real , Dim > > polygons( simplices.size() );
	for( unsigned int i=0 ; i<simplices.size() ; i++ )
	{
		polygons[i].resize( K+1 );
		for( unsigned int k=0 ; k<=K ; k++ ) polygons[i][k] = simplices[i][k];
		polygons[i].v = vf[i];
	}
	PLY::WritePolygons< Factory , PlyVFFace< Index , Real , Dim > >( fileName , Factory() , vertices , polygons , &properties[0] , static_cast< unsigned int >( properties.size() ) , file_type , nullptr );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , Dim > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , const std::vector< Point< Real , Dim > > &vf , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::NormalFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	if( simplices.size()!=vf.size() ) MK_THROW( "Number of simplices and vectors don't match: " , simplices.size() , " != " , vf.size() );
	if( vertices.size()!=normals.size() ) MK_THROW( "Number of vertices and normals don't match: " , vertices.size() , " != " , normals.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = normals[i];

	std::vector< GregTurk::PlyProperty > properties = PlyVFFace< Index , Real , Dim >::Properties();
	std::vector< PlyVFFace< Index , Real , Dim > > polygons( simplices.size() );
	for( unsigned int i=0 ; i<simplices.size() ; i++ )
	{
		polygons[i].resize( K+1 );
		for( unsigned int k=0 ; k<=K ; k++ ) polygons[i][k] = simplices[i][k];
		polygons[i].v = vf[i];
	}
	PLY::WritePolygons< Factory , PlyVFFace< Index , Real , Dim > >( fileName , Factory() , _vertices , polygons , &properties[0] , static_cast< unsigned int >( properties.size() ) , file_type , nullptr );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteColorVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , 3 > > &colors , const std::vector< SimplexIndex< K , Index > > &simplices , const std::vector< Point< Real , Dim > > &vf , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > >;
	using Vertex = typename Factory::DataType;

	if( simplices.size()!=vf.size() ) MK_THROW( "Number of simplices and vectors don't match: " , simplices.size() , " != " , vf.size() );
	if( vertices.size()!=colors.size() ) MK_THROW( "Number of vertices and colors don't match: " , vertices.size() , " != " , colors.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = colors[i]*255;

	std::vector< GregTurk::PlyProperty > properties = PlyVFFace< Index , Real , Dim >::Properties();
	std::vector< PlyVFFace< Index , Real , Dim > > polygons( simplices.size() );
	for( unsigned int i=0 ; i<simplices.size() ; i++ )
	{
		polygons[i].resize( K+1 );
		for( unsigned int k=0 ; k<=K ; k++ ) polygons[i][k] = simplices[i][k];
		polygons[i].v = vf[i];
	}
	PLY::WritePolygons< Factory , PlyVFFace< Index , Real , Dim > >( fileName , Factory() , _vertices , polygons , &properties[0] , static_cast< unsigned int >( properties.size() ) , file_type , nullptr );
}

template< typename Index , typename Real , unsigned int K , unsigned int Dim >
void WriteColorVectorFieldMesh( std::string fileName , const std::vector< Point< Real , Dim > > &vertices , const std::vector< Point< Real , 3 > > &colors , const std::vector< Point< Real , Dim > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , const std::vector< Point< Real , Dim > > &vf , int file_type )
{
	using Factory = DataFactory::Factory< Real , DataFactory::PositionFactory< Real , Dim > , DataFactory::RGBColorFactory< Real > , DataFactory::NormalFactory< Real , Dim > >;
	using Vertex = typename Factory::DataType;

	if( simplices.size()!=vf.size() ) MK_THROW( "Number of simplices and vectors don't match: " , simplices.size() , " != " , vf.size() );
	if( vertices.size()!=colors.size() ) MK_THROW( "Number of vertices and colors don't match: " , vertices.size() , " != " , colors.size() );
	if( vertices.size()!=normals.size() ) MK_THROW( "Number of vertices and normals don't match: " , vertices.size() , " != " , normals.size() );

	std::vector< Vertex > _vertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = colors[i]*255 , _vertices[i].template get<2>() = normals[i];

	std::vector< GregTurk::PlyProperty > properties = PlyVFFace< Index , Real , Dim >::Properties();
	std::vector< PlyVFFace< Index , Real , Dim > > polygons( simplices.size() );
	for( unsigned int i=0 ; i<simplices.size() ; i++ )
	{
		polygons[i].resize( K+1 );
		for( unsigned int k=0 ; k<=K ; k++ ) polygons[i][k] = simplices[i][k];
		polygons[i].v = vf[i];
	}
	PLY::WritePolygons< Factory , PlyVFFace< Index , Real , Dim > >( fileName , Factory() , _vertices , polygons , &properties[0] , static_cast< unsigned int >( properties.size() ) , file_type , nullptr );
}

//////////////////////////
// Helper Functionality //
//////////////////////////

template< typename Index , typename Real , unsigned int K >
std::vector< Point< Real , K+1 > > SimplexToVertexNormals( const std::vector< Point< Real , K+1 > > &vertices , const std::vector< SimplexIndex< K , Index > > &simplices )
{
	static const unsigned int Dim = K+1;
	std::vector< Point< Real , Dim > > normals( vertices.size() );
	for( unsigned int i=0 ; i<simplices.size() ; i++ )
	{
		Simplex< Real , Dim , K > s;
		for( unsigned int k=0 ; k<=K ; k++ ) s[k] = vertices[ simplices[i][k] ];
		Point< Real , Dim > n = s.normal();
		for( unsigned int k=0 ; k<=K ; k++ ) normals[ simplices[i][k] ] += n;
	}

	return normals;
}

///////////////
// PlyVFFace //
///////////////

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim >::PlyVFFace( void ) : vertices(nullptr) , nr_vertices(0) {}

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim >::~PlyVFFace( void ){ resize(0); }

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim >::PlyVFFace( const PlyVFFace & face )
{
	vertices = nullptr;
	(*this) = face;
}

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim > & PlyVFFace< Index , Real , Dim >::operator = ( const PlyVFFace& face )
{
	if( vertices ) free( vertices ) , vertices = nullptr;
	nr_vertices = face.nr_vertices;
	if( nr_vertices ) vertices = (Index*)malloc( sizeof(Index)*nr_vertices );
	else              vertices = nullptr;
	memcpy( vertices , face.vertices , sizeof(Index)*nr_vertices );
	v = face.v;
	return *this;
}

template< typename Index , typename Real , unsigned int Dim >
void PlyVFFace< Index , Real , Dim >::resize( unsigned int count )
{
	if( vertices ) free( vertices ) , vertices = nullptr;
	nr_vertices = 0;
	if( count ) vertices = (Index*)malloc( sizeof(Index)*count ) , nr_vertices = count;
}

template< typename Index , typename Real , unsigned int Dim >
unsigned int PlyVFFace< Index , Real , Dim >::size( void ) const { return nr_vertices; }

template< typename Index , typename Real , unsigned int Dim >
Index& PlyVFFace< Index , Real , Dim >::operator[] ( unsigned int idx ){ return vertices[idx]; }

template< typename Index , typename Real , unsigned int Dim >
const Index & PlyVFFace< Index , Real , Dim >::operator[] ( unsigned int idx ) const { return vertices[idx]; }

template< typename Index , typename Real , unsigned int Dim >
std::vector< GregTurk::PlyProperty > PlyVFFace< Index , Real , Dim >::Properties( void )
{
	std::vector< GregTurk::PlyProperty > properties( Dim+1 );
	properties[0] = GregTurk::PlyProperty( "vertex_indices" , PLY::Type< Index >() ,  PLY::Type< Index >() , (int)offsetof( PlyVFFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , nr_vertices ) );
	for( unsigned int d=0 ; d<Dim ; d++ ) properties[d+1] = GregTurk::PlyProperty( "vf_" + std::to_string(d) , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*d , 0 , 0 , 0 , 0 );
	return properties;
};

template< typename Index , typename Real , unsigned int Dim >
std::vector< GregTurk::PlyProperty > PlyVFFace< Index , Real , Dim >::ReadProperties( void )
{
	std::vector< GregTurk::PlyProperty > properties( Dim+1 );
	properties[0] = GregTurk::PlyProperty( "vertex_indices" , PLY::Type< Index >() ,  PLY::Type< Index >() , (int)offsetof( PlyVFFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , nr_vertices ) );
	for( unsigned int d=0 ; d<Dim ; d++ ) properties[d+1] = GregTurk::PlyProperty( "vf_" + std::to_string(d) , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*d , 0 , 0 , 0 , 0 );
	if( Dim>=1 ) properties.push_back( GregTurk::PlyProperty( "vx" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords[0] ) , 0 , 0 , 0 , 0 ) );
	if( Dim>=2 ) properties.push_back( GregTurk::PlyProperty( "vy" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords[1] ) , 0 , 0 , 0 , 0 ) );
	if( Dim>=3 ) properties.push_back( GregTurk::PlyProperty( "vz" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords[2] ) , 0 , 0 , 0 , 0 ) );
	if( Dim>=4 ) properties.push_back( GregTurk::PlyProperty( "vw" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords[3] ) , 0 , 0 , 0 , 0 ) );
	return properties;
};

/////////////////////
// PlyTexturedFace //
/////////////////////

template< typename Index , typename Real >
PlyTexturedFace< Index , Real >::PlyTexturedFace( void ){ vertices = nullptr , uv_coordinates = nullptr , nr_vertices = nr_uv_coordinates = 0; }

template< typename Index , typename Real >
PlyTexturedFace< Index , Real >::~PlyTexturedFace( void ){ resize(0); }

template< typename Index , typename Real >
PlyTexturedFace< Index , Real >::PlyTexturedFace( const PlyTexturedFace& face ) : vertices(nullptr) , uv_coordinates(nullptr) { *this = face; }

template< typename Index , typename Real >
PlyTexturedFace< Index , Real > & PlyTexturedFace< Index , Real >::operator = ( const PlyTexturedFace& face )
{
	if( vertices ) free( vertices ) , vertices = nullptr;
	if( uv_coordinates ) free( uv_coordinates ) , uv_coordinates = nullptr;
	nr_vertices = face.nr_vertices , nr_uv_coordinates = face.nr_uv_coordinates;
	if( nr_vertices ) vertices = (Index*)malloc( sizeof(Index)*nr_vertices );
	else              vertices = nullptr;
	if( nr_uv_coordinates ) uv_coordinates = (Real*)malloc( sizeof(Real)*nr_uv_coordinates );
	else                    uv_coordinates = nullptr;
	memcpy( vertices , face.vertices , sizeof(Index)*nr_vertices );
	memcpy( uv_coordinates , face.uv_coordinates , sizeof(Real)*nr_uv_coordinates );
	return *this;
}

template< typename Index , typename Real >
void PlyTexturedFace< Index , Real >::resize( unsigned int count )
{
	if( vertices ) free( vertices ) , vertices = nullptr;
	if( uv_coordinates ) free( uv_coordinates ) , uv_coordinates = nullptr;
	nr_vertices = nr_uv_coordinates = 0;
	if( count )
	{
		vertices = (Index*)malloc( sizeof(Index)*count ) , nr_vertices = count;
		uv_coordinates = (Real*)malloc( sizeof(Real)*2*count ) , nr_uv_coordinates = count*2;
	}
}

template< typename Index , typename Real >
Index& PlyTexturedFace< Index , Real >::operator[] ( unsigned int idx ){ return vertices[idx]; }

template< typename Index , typename Real >
const Index & PlyTexturedFace< Index , Real >::operator[] ( unsigned int idx ) const { return vertices[idx]; }

template< typename Index , typename Real >
const Point2D< Real > & PlyTexturedFace< Index , Real >::texture( unsigned int idx ) const { return Point2D< Real >( uv_coordinates[2*idx] , uv_coordinates[2*idx+1] ); }

template< typename Index , typename Real >
Point2D< Real > & PlyTexturedFace< Index , Real >::texture( unsigned int idx ){ return *( (Point2D< Real >*)(uv_coordinates+2*idx) ); }

template< typename Index , typename Real >
int PlyTexturedFace< Index , Real >::size( void ) const { return nr_vertices; }

template< typename Index , typename Real >
GregTurk::PlyProperty PlyTexturedFace< Index , Real >::ReadProperties[] =
{
	{ "vertex_indices" , PLY::Type< Index >() , PLY::Type< Index >() , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
	{ "texcoord" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyTexturedFace , uv_coordinates ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_uv_coordinates ) } ,
};

template< typename Index , typename Real >
GregTurk::PlyProperty PlyTexturedFace< Index , Real >::WriteProperties[] =
{
	{ "vertex_indices" , PLY::Type< Index >() , PLY::Type< Index >() , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_UCHAR , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
	{ "texcoord", PLY::Type< Real >(), PLY::Type< Real >(), (int)offsetof(PlyTexturedFace, uv_coordinates), 1, PLY_UCHAR, PLY_INT, (int)offsetof(PlyTexturedFace, nr_uv_coordinates) },
};

