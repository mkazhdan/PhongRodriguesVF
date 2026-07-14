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

#include <vector>

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"

#include "Include/PreProcessing.h"
#include "Include/PlyIO.h"
#include "Include/NormalFitter.h"
#include "Include/SimplexProcessing.h"
#include "Include/SimplicialMeshProcessing.h"

using namespace MishaK;
using namespace MishaK::SimplicialMesh;

static const unsigned int K = 2;
static const unsigned int Dim = K+1;

CmdLineParameterArray< std::string , 3 >
	In( "in" );

CmdLineParameter< std::string >
	Out( "out" );

CmdLineParameter< int >
	LoopNormalIterations( "loopIters" , LOOP_NORMAL_ITERS );

CmdLineParameter< unsigned int >
	QuadratureSamples( "qSamples" , DEFAULT_QUADRATURE );

CmdLineReadable
	CovariantDerivativeDifference( "covDiff" ) ,
	Verbose( "verbose" );


std::vector< CmdLineReadable * > params =
{
	&In ,
	&Out ,
	&QuadratureSamples ,
	&LoopNormalIterations ,
	&CovariantDerivativeDifference ,
	&Verbose ,
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh, vector field 1, vector field 2>\n" , In.name.c_str()  );
	printf( "\t[--%s <output vector field>]\n" , Out.name.c_str() );
	printf( "\t[--%s <quadrature samples in {%s}>=%d]\n" , QuadratureSamples.name.c_str() , QuadratureValuesString().c_str() , QuadratureSamples.value );
	printf( "\t[--%s <loop normal iterations>=%d]\n" , LoopNormalIterations.name.c_str() , LoopNormalIterations.value );
	printf( "\t[--%s]\n" , CovariantDerivativeDifference.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

template< unsigned int Quadrature >
std::vector< Point< double , Dim > > GetBracket( const EmbeddedPhongMesh< K > &mesh , const std::vector< Point< double , Dim > > & X , const std::vector< Point< double , Dim > > & Y )
{
	std::vector< Point< double , Dim > > Z( mesh.vertexNum() );

	Miscellany::Timer timer;

	Eigen::SparseMatrix< double > P = mesh.tangentProlongation();
	Eigen::SparseMatrix< double > Pt = P.transpose();

	Eigen::SparseMatrix< double > M = Pt * mesh.template massMatrix< Quadrature >() * P;
	if( Verbose.set ) std::cout << "Got system matrices: " << timer() << std::endl;

	LLtSolver solver( M );
	if( Verbose.set ) std::cout << "Factorized system matrices: " << timer() << std::endl;

	// Compute the bracket/difference of covariant derivatives of the two vector fields and integrate against the Phong-Rodrigues vector-field basis
	Eigen::VectorXd _Z;
	if( CovariantDerivativeDifference.set ) _Z = mesh.template massVector< Quadrature >( SimplicialMesh::PhongRodriguesCovariantDerivativeDifferenceField< K >( mesh , X , Y ) );
	else                                    _Z = mesh.template massVector< Quadrature >( SimplicialMesh::PhongRodriguesLieBracketField< K >( mesh , X , Y ) );
	if( Verbose.set )
	{
		if( CovariantDerivativeDifference.set ) std::cout << "Got difference of covariant derivatives: " << timer() << std::endl;
		else                                    std::cout << "Got Lie bracket: " << timer() << std::endl;
	}

	// Compute the least-squares best-fit within the space spanned by the Phong-Rodrigues vector-field basis
	_Z = P * solver.solve( Pt * _Z );
	if( Verbose.set ) std::cout << "Solved system: " << timer() << std::endl;

	for( size_t i=0 ; i<Z.size() ; i++ ) for( unsigned int k=0 ; k<=2 ; k++ ) Z[i][k] = _Z[i*3+k];

	return Z;
}

template< unsigned int Quadrature , unsigned int ... Quadratures >
std::vector< Point< double , Dim > > Execute( unsigned int quadrature , const EmbeddedPhongMesh< K > &mesh , const std::vector< Point< double , Dim > > & X , const std::vector< Point< double , Dim > > & Y )
{
	if( quadrature==Quadrature ) return GetBracket< Quadrature >( mesh , X , Y );
	else if constexpr( sizeof...(Quadratures) ) return Execute< Quadratures... >( quadrature , mesh , X , Y );
	else
	{
		MK_THROW( "Bad quadrature option: " , quadrature );
		return std::vector< Point< double , Dim > >();
	}
}

std::vector< Point< double , Dim > > ReadVectorField( std::string fileName )
{
	std::vector< Point< double , Dim > > vf;
	std::string ext = ToLower( GetFileExtension( fileName ) );
	if( ext=="ply" )
	{
		using Factory = DataFactory::VectorFieldFactory< double , Dim >;
		using Vertex = typename Factory::DataType;
		PLY::ReadVertices( fileName , Factory() , vf );
	}
	else
	{
		std::ifstream iStream( fileName );
		Point< double , Dim > p;
		while( iStream >> p[0] >> p[1] >> p[2] ) vf.push_back( p );
	}
	return vf;
}

void WriteVectorField( std::string fileName , const std::vector< Point< double , Dim > > & vertices , const std::vector< Point< double , Dim > > & normals , const std::vector< SimplexIndex< K > > & simplices , const std::vector< Point< double , Dim > > & vf , bool ascii )
{
	std::string ext = ToLower( GetFileExtension( fileName ) );
	if( ext=="ply" )
	{
		using Factory = DataFactory::Factory< double , DataFactory::PositionFactory< double , Dim > , DataFactory::NormalFactory< double , Dim > , DataFactory::VectorFieldFactory< double , Dim > >;
		using Vertex = typename Factory::DataType;
		std::vector< Vertex > _vertices( vertices.size() );
		for( size_t i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = normals[i] , _vertices[i].template get<2>() = vf[i];
		PLY::WriteTriangles( fileName , Factory() , _vertices , simplices , ascii ? PLY_ASCII : PLY_BINARY_NATIVE );
	}
	else
	{
		std::ofstream oStream( fileName );
		for( size_t i=0 ; i<vf.size() ; i++ ) oStream << vf[i][0] << " " << vf[i][1] << " " << vf[i][2] << std::endl;
	}
}

int main( int argc , char * argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	bool ascii = false;
	std::vector< Point< double , Dim > > vertices , normals;
	std::vector< SimplexIndex< K > > simplices;
	std::vector< Point3D< double > > X , Y , Z;

	// Get the mesh
	{
		ascii = PlyIO::ReadMesh( In.values[0] , vertices , normals , simplices ) == PLY_ASCII;

		// If per-vertex normals aren't provided, compute them from the incident triangles
		if( normals.size()!=vertices.size() )
		{
			if( Verbose.set ) std::cout << "Computing normals" << std::endl;
			typename NormalFitter< K >::Params params;
			if( LoopNormalIterations.value<0 ) params.subdivisionIterations = -LoopNormalIterations.value , params.useWarren = false;
			else                               params.subdivisionIterations =  LoopNormalIterations.value , params.useWarren = true;
			params.quadrature = QuadratureSamples.value;
			normals = NormalFitter< K >::Fit( vertices , simplices , params );
		}

		// Normalize the normals to have unit length
		for( unsigned int i=0 ; i<normals.size() ; i++ ) normals[i] = normals[i] /= Point< double , Dim >::Length( normals[i] );
	}

	// Get the vector fields
	{
		X = ReadVectorField( In.values[1] );
		Y = ReadVectorField( In.values[2] );

		if( X.size()!=vertices.size() ) MK_ERROR_OUT( "X vector-field size does not match vertex num: " , X.size() , " != " , vertices.size() );
		if( Y.size()!=vertices.size() ) MK_ERROR_OUT( "Y vector-field size does not match vertex num: " , Y.size() , " != " , vertices.size() );
	}

	// Subtract off the normal components from the vector-fields
	{
		for( unsigned int i=0 ; i<normals.size() ; i++ )
		{
			X[i] -= normals[i] * Point< double , Dim >::Dot( normals[i] , X[i] );
			Y[i] -= normals[i] * Point< double , Dim >::Dot( normals[i] , Y[i] );
		}
	}

	{
		EmbeddedPhongMesh< K > mesh( vertices , normals , simplices );

		try{ Z = Execute< SUPPORTED_QUADRATURE >( QuadratureSamples.value , mesh , X , Y ); }
		catch( Exception & ){ MK_ERROR_OUT( "Only quadrature samples in {" , QuadratureValuesString() , "} supported" ); }
	}

	if( Out.set ) WriteVectorField( Out.value , vertices , normals , simplices , Z , ascii );

	return EXIT_SUCCESS;
}
