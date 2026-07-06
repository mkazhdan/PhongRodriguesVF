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

#ifndef SIMPLICIAL_MESH_PROCESSING_CONCEPTS_INCLUDED
#define SIMPLICIAL_MESH_PROCESSING_CONCEPTS_INCLUDED


namespace MishaK
{
	namespace SimplicialMesh
	{
		template< typename Field , unsigned int K , typename T >
		concept HasMeshFunction = SimplexProcessing::HasArrayOfSimplexFunctions< Field , K , T >;

		template< typename Field , unsigned int K >
		concept HasMeshScaleFactorFunction = HasMeshFunction< Field , K , SimplexProcessing::ScaleFactor >;

		template< typename Field , unsigned int K >
		concept HasMeshTangentVectorFunction = HasMeshFunction< Field , K , SimplexProcessing::Differential< K , double > >;

		template< typename Field , unsigned int K , typename T >
		concept HasMeshLinearMapFunction = SimplexProcessing::HasArrayOfSimplexlinearMapFunctions< Field , K , T >;

		template< typename Field , unsigned int K , unsigned int N , typename T >
		concept HasMeshSystemFunction = SimplexProcessing::HasArrayOfSimplexSystemMatrixFunctions< Field , K , N , T >;

		template< typename Field , unsigned int K , typename T >
		concept HasMeshBilinearFormFunction = SimplexProcessing::HasArrayOfSimplexBilinearFormFunctions< Field , K , T >;

		template< typename Field , unsigned int K , unsigned int N , typename T >
		concept HasMeshSystemLinearMapOrBilinearFormFunction = HasMeshSystemFunction< Field , K , N , T > || HasMeshLinearMapFunction< Field , K , T > || HasMeshBilinearFormFunction< Field , K , T >;

		template< typename Field , unsigned int K , unsigned int N , typename T >
		concept HasMeshDSystemFunction = SimplexProcessing::HasArrayOfSimplexSystemMatrixFunctions< Field , K , N , SimplexProcessing::Differential< K , T > >;

		template< typename Field , unsigned int K , unsigned int N , typename T >
		concept HasMeshSystemAndDSystemFunction = SimplexProcessing::HasArrayOfSimplexSystemMatrixFunctions< Field , K , N , std::pair< T , SimplexProcessing::Differential< K , T > > >;

		template< typename IndexFunctor >
		concept HasIndexFunctor = requires( const IndexFunctor f , size_t sIdx , unsigned int eIdx ) { { f(sIdx,eIdx) } -> std::convertible_to< size_t >; };
	}
}
#endif // SIMPLICIAL_MESH_PROCESSING_CONCEPTS_INCLUDED
