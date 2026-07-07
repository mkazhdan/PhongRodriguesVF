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

///////////////////////
// SystemVectorField //
///////////////////////
template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexSystemVectorFunction< K , NumF , T > SystemVectors >
auto SystemVectorField( TestFunctions && Fs , SystemVectors && S )
{
	return [_Fs=std::forward< TestFunctions >(Fs),_S=std::forward< SystemVectors >(S)]( SimplexProcessing::Position< K > p )
		{
			T fs[NumF];
			for( unsigned int n=0 ; n<NumF ; n++ ) fs[n] = _Fs[n](p);
			return _S(p)( fs );
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexlinearMapFunction< K , T > LinearMaps >
auto SystemVectorField( Function && F , TestFunctions && Fs , LinearMaps && L )
{
	return [_F=std::forward< Function>(F),_Fs=std::forward< TestFunctions >(Fs),_L=std::forward< LinearMaps >(L)]( SimplexProcessing::Position< K > p )
		{
			Point< double , NumF > v;

			T f = _L(p)( _F(p) );
			for( unsigned int n=0 ; n<NumF ; n++ ) v[n] = DotProduct( f , _Fs[n](p) );
			return v;
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexBilinearFormFunction< K , T > BilinearForms >
	requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > )
auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B )
{
	return [_F=std::forward< Function >(F),_Fs=std::forward< TestFunctions >(Fs),_B=std::forward< BilinearForms >(B)]( SimplexProcessing::Position< K > p )
		{
			Point< double , NumF > v;

			T f = _F(p);
			const auto & b = _B(p);
			for( unsigned int n=0 ; n<NumF ; n++ ) v[n] = b( f , _Fs[n](p) );
			return v;
		};
}

///////////////////////
// SystemMatrixField //
///////////////////////
template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexSystemMatrixFunction< K , NumF , T > SystemMatrices >
	requires( !HasSimplexlinearMapFunction< SystemMatrices , K , T > && !HasSimplexBilinearFormFunction< SystemMatrices , K , T > )
auto SystemMatrixField( TestFunctions && Fs , SystemMatrices && S )
{
	return [_Fs=std::forward< TestFunctions >(Fs),_S=std::forward< SystemMatrices >(S)]( SimplexProcessing::Position< K > p )
		{
			T fs[NumF];
			for( unsigned int n=0 ; n<NumF ; n++ ) fs[n] = _Fs[n](p);
			return _S(p)( fs );
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexlinearMapFunction< K , T > LinearMaps >
auto SystemMatrixField( TestFunctions && Fs , LinearMaps && L )
{
	return [_Fs=std::forward< TestFunctions >(Fs),_L=std::forward< LinearMaps >(L)]( SimplexProcessing::Position< K > p )
		{
			SquareMatrix< double , NumF > m;
			T fs[NumF] , lfs[NumF];

			const auto & l = _L(p);
			for( unsigned int n=0 ; n<NumF ; n++ )
			{
				fs[n] = _Fs[n](p);
				lfs[n] = l( fs[n] );
			}
			for( unsigned int n1=0 ; n1<NumF ; n1++ ) for( unsigned int n2=0 ; n2<NumF; n2++ ) m(n1,n2) = DotProduct( fs[n1] , lfs[n2] );
			return m;
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexBilinearFormFunction< K , T > BilinearForms >
	requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > )
auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B )
{
	return [_Fs=std::forward< TestFunctions >(Fs),_B=std::forward< BilinearForms >(B)]( SimplexProcessing::Position< K > p )
		{
			SquareMatrix< double , NumF > m;
			T fs[NumF];

			const auto & b = _B(p);
			for( unsigned int n=0 ; n<NumF ; n++ ) fs[n] = _Fs[n](p);
			for( unsigned int n1=0 ; n1<NumF ; n1++ ) for( unsigned int n2=0 ; n2<NumF; n2++ ) m(n1,n2) = b( fs[n1] , fs[n2] );
			return m;
		};
}

//////////////////
// MCIntegrator //
//////////////////

template< unsigned int K , unsigned int QuadratureSamples >
template< typename T , HasSimplexFunction< K , T > Function >
T MCIntegrator< K , QuadratureSamples >::Integral( Function && F )
{
	return SimplexProcessing::SampleSum< T >( ArrayWrapper( [&]( size_t q ){ return F( Position( static_cast< unsigned int >(q) ) ) * Measure( static_cast< unsigned int >(q) ); } ) , QuadratureSamples );
}
