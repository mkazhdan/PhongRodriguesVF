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
#include <set>

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/PlyData.h"

#include "Include/PreProcessing.h"
#include "Include/SimplicialMeshProcessing.h"
#include "Include/PlyIO.h"
#include "Include/NormalFitter.h"


using namespace MishaK;
using namespace MishaK::SimplicialMesh;

static const unsigned int K = 2;
static const unsigned int Dim = K+1;

enum EnergyType
{
	Connection ,
	Hodge ,
	Killing
};

const std::vector< std::string > EnergyTypeNames = { "Connection" , "Hodge" , "Killing" };

CmdLineParameterArray< std::string , 2 >
	In( "in" );

CmdLineParameter< std::string >
	Out( "out" );

CmdLineParameter< int >
	LoopNormalIterations( "loopIters" , LOOP_NORMAL_ITERS );

CmdLineParameter< unsigned int >
	Energy( "energy" , EnergyType::Connection ) ,
	QuadratureSamples( "qSamples" , DEFAULT_QUADRATURE );

CmdLineReadable
	Verbose( "verbose" );

std::vector< CmdLineReadable * > params =
{
	&In ,
	&Out ,
	&Verbose ,
	&QuadratureSamples ,
	&Energy ,
	&LoopNormalIterations ,
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh, input constraints>\n" , In.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str()  );
	printf( "\t[--%s <quadrature samples in {%s}>=%d]\n" , QuadratureSamples.name.c_str() , QuadratureValuesString().c_str() , QuadratureSamples.value );
	printf( "\t[--%s <energy type>=%d]\n" , Energy.name.c_str() , Energy.value );
	for( unsigned int i=0 ; i<EnergyTypeNames.size() ; i++ ) printf( "\t\t%d] %s\n" , i , EnergyTypeNames[i].c_str() );
	printf( "\t[--%s <loop normal iterations>=%d]\n" , LoopNormalIterations.name.c_str() , LoopNormalIterations.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

template< unsigned int Quadrature >
std::vector< Point< double , Dim > > FitVF
(
	EmbeddedPhongMesh< K > & mesh ,
	const std::vector< std::pair< size_t , Point< double , Dim > > > & constraints
)
{
	Miscellany::Timer timer;

	std::vector< Point< double , Dim > > outVF( mesh.vertexNum() );
	size_t dofs = mesh.vertexNum() - constraints.size();

	// The prolongation matrix mapping from degrees of freedom to the tangent vector coefficients
	Eigen::SparseMatrix< double > P1( mesh.vertexNum() * 2 , dofs * 2 );
	{
		std::vector< Eigen::Triplet< double > > triplets( dofs*2 );

		std::set< size_t > constrainedIndices;
		for( unsigned int i=0 ; i<constraints.size() ; i++ ) constrainedIndices.insert( constraints[i].first );

		size_t idx = 0;
		for( unsigned int i=0 ; i<mesh.vertexNum() ; i++ ) if( constrainedIndices.find( i )==constrainedIndices.end() ) 
		{
			triplets[2*idx+0] = Eigen::Triplet< double >( 2*i+0 , static_cast< unsigned int >( 2*idx+0 ) , 1. );
			triplets[2*idx+1] = Eigen::Triplet< double >( 2*i+1 , static_cast< unsigned int >( 2*idx+1 ) , 1. );
			idx++;
		}

		P1.setFromTriplets( triplets.begin() , triplets.end() );
	}

	// The prolongation matrix from framed tangent vector coefficients to extrinsic tangent vector coeffidients
	Eigen::SparseMatrix< double > P2 = mesh.tangentProlongation();

	// The stiffness matrix
	Eigen::SparseMatrix< double > S;
	switch( Energy.value )
	{
	case EnergyType::Connection: S = mesh.template stiffness< Quadrature >                                () ; break;
	case EnergyType::Hodge:      S = mesh.template stiffness< Quadrature , CovariantComponent::Hodge >    () ; break;
	case EnergyType::Killing:    S = mesh.template stiffness< Quadrature , CovariantComponent::Symmetric >() ; break;
	default: MK_ERROR_OUT( "Unrecognized energy type: " , Energy.value );
	}

	Eigen::SparseMatrix< double > P1t = P1.transpose() , P2t = P2.transpose();
	Eigen::VectorXd C( mesh.vertexNum()*Dim );

	C.setZero();
	for( unsigned int i=0 ; i<constraints.size() ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) C( constraints[i].first * Dim + d ) = constraints[i].second[d];
	if( Verbose.set ) std::cout << "Got system: " << timer() << std::endl;

	// Solve for the minimzer of:
	//	E(X) = ( P2*P1*X + C )^t * S * ( P2*P1*X + C )
	//       = X^T * P1^t * P2^t * S * P2 * P1 * X + 2 * X^T * P1^t * PT^2 * S * C + ...
	// => P1^t * P2^t * S * P2 * P1 * X = - P1^t * PT^2 * S * C

	timer.reset();
	LLtSolver solver( P1t * P2t * S * P2 * P1 );
	if( solver.info()!=Eigen::Success ) MK_ERROR_OUT( "Failed to factorize system" );
	if( Verbose.set ) std::cout << "Factored system: " << timer() << std::endl;

	timer.reset();
	Eigen::VectorXd x = C + P2 * P1 * solver.solve( - P1t * P2t * S * C );
	for( unsigned int i=0 ; i<mesh.vertexNum() ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) outVF[i][d] = x[i*Dim+d];
	if( Verbose.set ) std::cout << "Solved system: " << timer() << std::endl;

	return outVF;
}

template< unsigned int Quadrature , unsigned int ... Quadratures >
std::vector< Point< double , Dim > > FitVF
(
	unsigned int quadrature ,
	EmbeddedPhongMesh< K > & mesh ,
	const std::vector< std::pair< size_t , Point< double , Dim > > > & constraints
)
{
	if( quadrature==Quadrature ) return FitVF< Quadrature >( mesh , constraints );
	else if constexpr( sizeof...(Quadratures) ) return FitVF< Quadratures... >( quadrature , mesh , constraints );
	else
	{
		MK_THROW( "Bad quadrature option: " , quadrature );
		return std::vector< Point< double , Dim > >();
	}
}

std::vector< std::pair< size_t , Point< double , Dim > > > ReadConstraints( std::string fileName )
{
	std::vector< std::pair< size_t , Point< double , Dim > > > constraints;
	std::ifstream iStream( fileName );
	size_t idx;
	Point< double , Dim > p;
	while( iStream >> idx >> p[0] >> p[1] >> p[2] ) constraints.emplace_back( idx , p );
	return constraints;
}

void WriteVectorField( std::string fileName , const std::vector< Point< double , Dim > > & vertices , const std::vector< Point< double , Dim > > & normals , const std::vector< SimplexIndex< K > >  & triangles , const std::vector< Point< double , Dim > > & vf , bool ascii )
{
	std::string ext = ToLower( GetFileExtension( fileName ) );
	if( ext=="ply" )
	{
		using Factory = DataFactory::Factory< double , DataFactory::PositionFactory< double , Dim > , DataFactory::NormalFactory< double , Dim > , DataFactory::VectorFieldFactory< double , Dim > >;
		using Vertex = typename Factory::DataType;
		std::vector< Vertex > _vertices( vertices.size() );
		for( size_t i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i] , _vertices[i].template get<1>() = normals[i] , _vertices[i].template get<2>() = vf[i];
		PLY::WriteTriangles( fileName , Factory() , _vertices , triangles , ascii ? PLY_ASCII : PLY_BINARY_NATIVE );
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
	std::vector< Point< double , Dim > > vertices , normals , vf;
	std::vector< SimplexIndex< K > > simplices;
	std::vector< std::pair< size_t , Point< double , Dim > > > constraints;

	// Read the mesh
	{
		ascii = PlyIO::ReadMesh( In.values[0] , vertices , normals , simplices ) == PLY_ASCII;

		// If per-vertex normals aren't provided, compute them from the incident triangles
		if( normals.size()!=vertices.size() )
		{
			typename NormalFitter< K >::Params params;
			if( LoopNormalIterations.value<0 ) params.subdivisionIterations = -LoopNormalIterations.value , params.useWarren = false;
			else                               params.subdivisionIterations =  LoopNormalIterations.value , params.useWarren = true;
			params.quadrature = 3;
			normals = NormalFitter< K >::Fit( vertices , simplices , params );
		}

		// Normalize the normals to have unit length
		for( unsigned int i=0 ; i<normals.size() ; i++ ) normals[i] = normals[i] /= Point< double , Dim >::Length( normals[i] );
	}

	// Read the constraints
	{
		constraints = ReadConstraints( In.values[1] );
		if( !constraints.size() ) MK_ERROR_OUT( "Need at least one constraint vertex" );

		// Subtract off the normal components from the constraints (and check that indices are within bounds)
		for( unsigned int i=0 ; i<constraints.size() ; i++ )
			if( constraints[i].first>=normals.size() ) MK_ERROR_OUT( "Constrain index " , i , " out of bounds: " , constraints[i].first , " >= " , normals.size() );
			else constraints[i].second -= normals[i] * Point< double , Dim >::Dot( normals[i] , constraints[i].second );
	}

	EmbeddedPhongMesh< K > mesh( vertices , normals , simplices );

	if( Verbose.set )
	{
		std::cout << "Vertices / Simplices: " << vertices.size() << " / " << simplices.size() << std::endl;
		std::cout << "Constraints: " << constraints.size() << std::endl;
		for( unsigned int i=0 ; i<constraints.size() ; i++ ) std::cout << "\t" << constraints[i].first << " : " << constraints[i].second << std::endl;
	}

	try{ vf = FitVF< SUPPORTED_QUADRATURE >( QuadratureSamples.value , mesh , constraints ); }
	catch( Exception & ){ MK_ERROR_OUT( "Only quadrature samples in {" , QuadratureValuesString() , "} supported" ); }

	if( Out.set ) WriteVectorField( Out.value , vertices , normals , simplices , vf , ascii );

	return EXIT_SUCCESS;
}
