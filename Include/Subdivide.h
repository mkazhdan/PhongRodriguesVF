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

#ifndef SUBDIVIDE_INCLUDED
#define SUBDIVIDE_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include "Misha/Geometry.h"
#include "Misha/MultiThreading.h"

#define USE_ADJACENCY

namespace MishaK
{
	namespace Subdivide
	{
		//////////////////
		// Declarations //
		//////////////////
		enum struct SubdivisionType
		{
			PLANAR ,		// One-to-four in-plane subdivision
			LOOP_LOOP ,		// With Loop definition of beta
			LOOP_WARREN ,	// With Warren definition of beta
			BUTTERFLY		// With weight w=1/16
		};

#ifdef USE_ADJACENCY
		struct Adjacency
		{
			static const unsigned int K=2;
			using Edge = std::pair< size_t , size_t >;

			static size_t CornerIndex( size_t t , unsigned int k ){ return t*(K+1)+( k%(K+1) ); };
			static void FactorCornerIndex( size_t c , size_t & t , unsigned int & k ){ t = c/(K+1) , k = c%(K+1); };

			Adjacency( const std::vector< SimplexIndex< K > > & triangles , size_t vNum );

			size_t vNum( void ) const { return _vNum; }
			size_t eNum( void ) const { return _edges.size(); }
			size_t tNum( void ) const { return _triangles.size(); }

			Edge edge( size_t e ) const { return _edges[e]; }
			Edge dualEdgeVertices( size_t e ) const { return Edge( cornerToVertex( _dualEdgeCorners[e].first ) , cornerToVertex( _dualEdgeCorners[e].second ) ); }
			Edge dualEdgeCorners( size_t e ) const { return Edge( _dualEdgeCorners[e].first , _dualEdgeCorners[e].second ); }
			SimplexIndex< 2 > triangle( size_t t ) const { return _triangles[t]; }

			std::vector< size_t > oneRingCorners( size_t v ) const;
			std::vector< size_t > oneRingVertices( size_t v ) const;

			size_t cornerToEdge( size_t c ) const { return _cornerToEdge[c]; }
			size_t oppositeCorner( size_t c ) const;
			size_t cornerToVertex( size_t c ) const;

		protected:
			size_t _vNum;
			const std::vector< SimplexIndex< K > > & _triangles;
			std::vector< Edge > _edges , _dualEdgeCorners;
			std::vector< size_t > _vertexToCorner;
			std::vector< size_t > _cornerToEdge;
		};
#endif // USE_ADJACENCY

		struct Subdivider
		{
			static const unsigned int K = 2;

			Subdivider( const std::vector< SimplexIndex< K > > & simplices , size_t vNum );

			std::vector< SimplexIndex< K > > operator()( void ) const;

			template< typename T >
			std::vector< T > operator()( const std::vector< T > & vertexData , SubdivisionType type ) const;

			template< typename T >
			std::vector< T > operator()( const std::vector< T > & cornerData ) const;

			static Eigen::MatrixXd LoopStencil( unsigned int valence , bool useWarren );
		protected:
#ifdef USE_ADJACENCY
			Adjacency _adjacency;
#else // !USE_ADJACENCY
			using Edge = std::pair< unsigned int , unsigned int >;
			struct _EdgeData
			{
				_EdgeData( unsigned int idx=-1 , unsigned int c1=-1 , unsigned int c2=-1 ) : idx(idx) , c1(c1) , c2(c2){}
				unsigned int idx , c1 , c2;
			};

			size_t _vNum;
			const std::vector< SimplexIndex< K > > & _triangles;
			std::map< Edge , _EdgeData > _edges;
#endif // USE_ADJACENCY
			std::vector< unsigned int > _valence;
		};

		template< unsigned int Dim >
		void Subdivide( std::vector< Point< double , Dim > > & vertices , std::vector< SimplexIndex< 2 > > & simplices , SubdivisionType type=SubdivisionType::LOOP_LOOP );

		template< unsigned int Dim >
		void Subdivide( std::vector< Point< double , Dim > > & vertices , std::vector< Point< double , Dim > > & normals , std::vector< SimplexIndex< 2 > > & simplices , SubdivisionType type=SubdivisionType::LOOP_LOOP );
#include "Subdivide.inl"
	}
}
#endif // SUBDIVIDE_INCLUDED