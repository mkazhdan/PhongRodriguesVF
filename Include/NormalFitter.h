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

#ifndef NORMAL_FITTER_INCLUDED
#define NORMAL_FITTER_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <map>

#include "Misha/Geometry.h"

#include "Include/Subdivide.h"
#include "Include/SimplicialMeshProcessing.h"

namespace MishaK
{
	namespace SimplicialMesh
	{
		template< unsigned int K >
		struct NormalFitter
		{
			static const unsigned int Dim = K+1;
			struct Params
			{
				unsigned int subdivisionIterations , quadrature;
				bool useWarren , normalizeScale , scaleByAverageEdgeLength;
				double diffusionTime , spectralEps;

				Params( void ) : subdivisionIterations(0) , useWarren(false) , diffusionTime(0.) , normalizeScale(true) , spectralEps(1e-8) , quadrature(1) , scaleByAverageEdgeLength(true) {}
			};

			static std::vector< Point< double , K+1 > > Fit( std::vector< Point< double , Dim > > & vertices , std::vector< SimplexIndex< K > > & simplices , Params params );

		protected:
			static Eigen::MatrixXd _NormalLimitStencil( Eigen::MatrixXd stencil , double eps );

			static std::vector< Point< double , K+1 > > _NormalsFromSimplices( const std::vector< Point< double , Dim > > & vertices , const std::vector< SimplexIndex< K > > & simplices , unsigned int iters , bool useWarren );

			static std::vector< Point< double , K+1 > > _LimitNormalsFromSimplices( const std::vector< Point< double , Dim > > & vertices , const std::vector< SimplexIndex< K > > & simplices , bool useWarren , double eps );

			template< unsigned int Quadrature >
			static void _SmoothNormals( EmbeddedMesh< K , Dim > mesh , std::vector< Point< double , Dim > > & normals , double diffusionTime );

			template< unsigned int Quadrature , unsigned int ... Quadratures >
			static void _SmoothNormals( unsigned int quadrature , EmbeddedMesh< K , Dim > mesh , std::vector< Point< double , Dim > > & normals , double diffusionTime );

		};
#include "NormalFitter.inl"
	}
}
#endif // LOOP_NORMALS_INCLUDED