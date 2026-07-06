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

#ifndef PRE_PROCESSING_INCLUDED
#define PRE_PROCESSING_INCLUDED

#if _WIN32 || _WIN64
#define NOMINMAX
#endif // _WIN32 || _WIN64

#define NEW_CODE							// New/experimental code
#define FAST_COMPILE						// Remove some of the templating options


#define USE_MY_STEIN						// Use self-implemented Stein system

////////////////////////
// Default Parameters //
////////////////////////

#define LOOP_NORMAL_ITERS 10

#ifdef FAST_COMPILE
#define SUPPORTED_QUADRATURE 1 , 3 , 6 , 12 , 24 , 32
#else // !FAST_COMPILE
#define SUPPORTED_QUADRATURE 1 , 3 , 4 , 6 , 7 , 12 , 13 , 24 , 27 , 32
#endif // FAST_COMPILE

#define DEFAULT_QUADRATURE 3


////////////////////////////////
////////////////////////////////
////////////////////////////////

#include <sstream>
template< unsigned int Quadrature , unsigned int ... Quadratures >
void _SetQuadratureValuesStringStream( std::stringstream & sStream )
{
	if constexpr( sizeof ... ( Quadratures ) )
	{
		sStream << Quadrature << ",";
		_SetQuadratureValuesStringStream< Quadratures ... >( sStream );
	}
	else sStream << Quadrature;
}

inline std::string QuadratureValuesString( void )
{
	std::stringstream sStream;
	_SetQuadratureValuesStringStream< SUPPORTED_QUADRATURE >( sStream );
	return sStream.str();
}

#include <Eigen/Sparse>
using LLtSolver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
using LDLtSolver = Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >;


#endif // PRE_PROCESSING_INCLUDED
