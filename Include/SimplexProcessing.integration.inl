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

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexlinearMapFunction< K , T > BilinearForms >
auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B )
{
	return [F,Fs,B]( SimplexProcessing::Position< K > p )
		{
			Point< double , NumF > v;

			T f = B(p)( F(p) );
			for( unsigned int n=0 ; n<NumF ; n++ ) v[n] = DotProduct( f , Fs[n](p) );
			return v;
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexBilinearFormFunction< K , T > BilinearForms >
	requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > )
auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B )
{
	return [F,Fs,B]( SimplexProcessing::Position< K > p )
		{
			Point< double , NumF > v;

			T f = F(p);
			const auto & b = B(p);
			for( unsigned int n=0 ; n<NumF ; n++ ) v[n] = b( f , Fs[n](p) );
			return v;
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexSystemVectorFunction< K , NumF , T > BilinearForms >
	requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > && !HasSimplexBilinearFormFunction< BilinearForms , K , T > )
auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B )
{
	return [Fs,B]( SimplexProcessing::Position< K > p )
		{
			T fs[NumF];
			for( unsigned int n=0 ; n<NumF ; n++ ) fs[n] = Fs[n](p);
			return B(p)( fs );
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexlinearMapFunction< K , T > BilinearForms >
auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B )
{
	return [Fs,B]( SimplexProcessing::Position< K > p )
		{
			SquareMatrix< double , NumF > m;
			T fs[NumF] , bfs[NumF];

			const auto & b = B(p);
			for( unsigned int n=0 ; n<NumF ; n++ )
			{
				fs[n] = Fs[n](p);
				bfs[n] = b( fs[n] );
			}
			for( unsigned int n1=0 ; n1<NumF ; n1++ ) for( unsigned int n2=0 ; n2<NumF; n2++ ) m(n1,n2) = DotProduct( fs[n1] , bfs[n2] );
			return m;
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexBilinearFormFunction< K , T > BilinearForms >
	requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > )
auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B )
{
	return [Fs,B]( SimplexProcessing::Position< K > p )
		{
			SquareMatrix< double , NumF > m;
			T fs[NumF];

			const auto & b = B(p);
			for( unsigned int n=0 ; n<NumF ; n++ ) fs[n] = Fs[n](p);
			for( unsigned int n1=0 ; n1<NumF ; n1++ ) for( unsigned int n2=0 ; n2<NumF; n2++ ) m(n1,n2) = b( fs[n1] , fs[n2] );
			return m;
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasArrayOfSimplexFunctions< K , T > TestFunctions , HasSimplexSystemMatrixFunction< K , NumF , T > BilinearForms >
	requires( !HasSimplexlinearMapFunction< BilinearForms , K , T > && !HasSimplexBilinearFormFunction< BilinearForms , K , T > )
auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B )
{
	return [Fs,B]( SimplexProcessing::Position< K > p )
		{
			T fs[NumF];
			for( unsigned int n=0 ; n<NumF ; n++ ) fs[n] = Fs[n](p);
			return B(p)( fs );
		};
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , ScaleFactor > ScaleFactors , HasSimplexFunction< K , T > Function , HasArrayOfSimplexFunctions< K , T > TestFunctions , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
auto ScaledSystemVectorField( ScaleFactors && SF , Function && F , TestFunctions && Fs , BilinearForms && B )
{
	auto VF = SystemVectorField< K , NumF , T >( std::forward< Function >( F ) , std::forward< TestFunctions >( Fs ) , std::forward< BilinearForms >( B ) );
	return [VF,SF]( SimplexProcessing::Position< K > p ){ return VF(p) * SF(p); };
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , ScaleFactor > ScaleFactors , HasArrayOfSimplexFunctions< K , T > TestFunctions , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
auto ScaledSystemMatrixField( ScaleFactors && SF , TestFunctions && Fs , BilinearForms && B )
{
	auto MF = SystemMatrixField< K , NumF , T >( std::forward< TestFunctions >( Fs ) , std::forward< BilinearForms >( B ) );
	return [MF,SF]( SimplexProcessing::Position< K > p ){ return MF(p) * SF(p); };
}

template< unsigned int K , unsigned int NumF , typename T , HasSimplexFunction< K , ScaleFactor > ScaleFactors , HasSimplexFunction< K , T > Field >
auto ScaledField( ScaleFactors && SF , Field && F )
{
	return [SF,F]( SimplexProcessing::Position< K > p ){ return SF(p) * F(p); };
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
