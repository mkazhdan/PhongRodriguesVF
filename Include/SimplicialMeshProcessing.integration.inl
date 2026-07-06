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

template< unsigned int K , unsigned int NumF , typename T , typename Function /* = Field< K , T > */ , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
auto SystemVectorField( Function && F , TestFunctions && Fs , BilinearForms && B )
{
	return SimplexProcessing::ArrayWrapper( [F,Fs,B]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::SystemVectorField< K , NumF , T >( F[sIdx] , Fs[sIdx] , B[sIdx] ); } );
}

template< unsigned int K , unsigned int NumF , typename T , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
auto SystemMatrixField( TestFunctions && Fs , BilinearForms && B )
{
	return SimplexProcessing::ArrayWrapper( [Fs,B]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::SystemMatrixField< K , NumF , T >( Fs[sIdx] , B[sIdx] ); } );
}

template< unsigned int K , unsigned int NumF , typename T , typename ScaleFactors /* = Samples< Field< K , ScaleFactor > > */ , typename Function /* = Field< K , T > */ , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
auto ScaledSystemVectorField( ScaleFactors && SF , Function && F , TestFunctions && Fs , BilinearForms && B )
{
	return SimplexProcessing::ArrayWrapper( [SF,F,Fs,B]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::ScaledSystemVectorField< K , NumF , T >( SF[sIdx] , F[sIdx] , Fs[sIdx] , B[sIdx] ); } );
}

template< unsigned int K , unsigned int NumF , typename T , typename ScaleFactors /* = Samples< Field< K , ScaleFactor > > */ , typename TestFunctions /* = Samples< Fields< K , T > > */ , typename BilinearForms /* = Field< K , BilinearForm< T > */ >
auto ScaledSystemMatrixField( ScaleFactors && SF , TestFunctions && Fs , BilinearForms && B )
{
	return SimplexProcessing::ArrayWrapper( [SF,Fs,B]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::ScaledSystemMatrixField< K , NumF , T >( SF[sIdx] , Fs[sIdx] , B[sIdx] ); } );
}

template< unsigned int K , unsigned int NumF , typename T , typename ScaleFactors /* = Samples< Field< K , ScaleFactor > > */ , typename Field /* = Samples< K , T > */ >
auto ScaledField( ScaleFactors && SF , Field && F )
{
	return SimplexProcessing::ArrayWrapper( [SF,F]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::ScaledField< K , NumF , T >( SF[sIdx] , F[sIdx] ); } );
}

//////////////////
// MCIntegrator //
//////////////////

template< unsigned int K , unsigned int QuadratureSamples >
template< typename T , HasMeshFunction< K , T > Function  >
T MCIntegrator< K , QuadratureSamples >::Integral( size_t numS , Function && F )
{
	return SystemAssembler::Integral< T >
		(
			numS ,
			SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< T >( F[sIdx] ); } )
		);
}

template< unsigned int K , unsigned int QuadratureSamples >
template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , HasMeshFunction< K , Point< double , NumElementsPerSimplex > > VectorFunctor >
Eigen::VectorXd MCIntegrator< K , QuadratureSamples >::Vector( size_t numF , size_t numS , Index && Idx , VectorFunctor && VF )
{
	return SystemAssembler::Vector< NumElementsPerSimplex >
		(
			numF ,
			numS ,
			std::forward< Index >( Idx ) ,
			SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< Point< double , NumElementsPerSimplex > >( VF[sIdx] ); } )
		);
}

template< unsigned int K , unsigned int QuadratureSamples >
template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
Eigen::SparseMatrix< double > MCIntegrator< K , QuadratureSamples >::Matrix( size_t numF , size_t numS , Index && Idx , MatrixFunctor && MF )
{
	return SystemAssembler::Matrix< NumElementsPerSimplex >
		(
			numF ,
			numS ,
			std::forward< Index >( Idx ) ,
			SimplexProcessing::ArrayWrapper( [MF]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] ); } )
		);
}

template< unsigned int K , unsigned int QuadratureSamples >
template< unsigned int NumElementsPerSimplex , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
void MCIntegrator< K , QuadratureSamples >::SetMatrixEntries( EigenMatrixEntries< NumElementsPerSimplex > &eme , size_t numS , MatrixFunctor && MF )
{
	return SystemAssembler::SetMatrixEntries< NumElementsPerSimplex >
		(
			eme ,
			numS ,
			SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] ); } )
		);
}

/////////////////////
// SystemAssembler //
/////////////////////
template< typename T , SimplexProcessing::HasArray< T > Integrand >
T SystemAssembler::Integral( size_t numS , Integrand && integrand )
{
	std::vector< T > summands( numS );
	ThreadPool::ParallelFor( 0 , numS , [&]( size_t sIdx ){ summands[sIdx] = integrand[ sIdx ]; } );

	T I = {};
	for( size_t sIdx=0 ; sIdx<numS ; sIdx++ ) I += summands[ sIdx ];
	return I;
}

template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , SimplexProcessing::HasArray< Point< double , NumElementsPerSimplex > > VectorFunctor >
Eigen::VectorXd SystemAssembler::Vector( size_t numF , size_t numS , Index && Idx , VectorFunctor && VF )
{
	std::vector< double > dEntries( numS * NumElementsPerSimplex );

	ThreadPool::ParallelFor
	(
		0 , numS ,
		[&]( size_t sIdx )
		{
			Point< double , NumElementsPerSimplex > D = VF[ sIdx ];
			for( unsigned int e=0 ; e<NumElementsPerSimplex ; e++ ) dEntries[ sIdx*NumElementsPerSimplex + e ] = D[e];
		}
	);

	Eigen::VectorXd D( numF );
	for( size_t i=0 ; i<numF ; i++ ) D[i] = 0;
	for( size_t sIdx=0 ; sIdx<numS ; sIdx++ ) for( unsigned int e=0 ; e<NumElementsPerSimplex ; e++ ) D[ Idx(sIdx,e) ] += dEntries[ sIdx*NumElementsPerSimplex + e ];

	return D;
}

template< unsigned int NumElementsPerSimplex , HasIndexFunctor Index , SimplexProcessing::HasArray< SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
Eigen::SparseMatrix< double > SystemAssembler::Matrix( size_t numF , size_t numS , Index && Idx , MatrixFunctor && MF )
{
	std::vector< Eigen::Triplet< double > > triplets( numS * NumElementsPerSimplex * NumElementsPerSimplex );

	ThreadPool::ParallelFor
	(
		0 , numS ,
		[&]( size_t sIdx )
		{
			SquareMatrix< double , NumElementsPerSimplex > M = MF[ sIdx ];
			for( unsigned int e1=0 ; e1<NumElementsPerSimplex ; e1++ ) for( unsigned int e2=0 ; e2<NumElementsPerSimplex ; e2++ )
				triplets[ sIdx * NumElementsPerSimplex * NumElementsPerSimplex + e1 * NumElementsPerSimplex + e2 ] = Eigen::Triplet< double >( Idx(sIdx,e1) , Idx(sIdx,e2) , M(e1,e2) );
		}
	);
	Eigen::SparseMatrix< double > M( numF , numF );
	M.setFromTriplets( triplets.begin() , triplets.end() );

	return M;
}

template< unsigned int NumElementsPerSimplex , SimplexProcessing::HasArray< SquareMatrix< double , NumElementsPerSimplex > > MatrixFunctor >
void SystemAssembler::SetMatrixEntries( EigenMatrixEntries< NumElementsPerSimplex > &eme , size_t numS , MatrixFunctor && MF )
{
	eme.clear();

	ThreadPool::ParallelFor
	(
		0 , numS ,
		[&]( size_t sIdx )
		{
			size_t idx = sIdx * NumElementsPerSimplex * NumElementsPerSimplex;
			SquareMatrix< double , NumElementsPerSimplex > M = MF[ sIdx ];
			for( unsigned int e1=0 ; e1<NumElementsPerSimplex ; e1++ ) for( unsigned int e2=0 ; e2<NumElementsPerSimplex ; e2++ , idx++ )
				MishaK::AddAtomic( eme.matrixEntry( idx ) , M(e1,e2) );
		}
	);
}

////////////////////////
// EigenMatrixEntries //
////////////////////////
template< unsigned int NumElementsPerSimplex >
template< HasIndexFunctor Index >
EigenMatrixEntries< NumElementsPerSimplex >::EigenMatrixEntries( size_t numF , size_t numS , Index && Idx ) :
	_M( new Eigen::SparseMatrix< double >( numF , numF ) ) , _matrixEntries( numS * NumElementsPerSimplex * NumElementsPerSimplex )
{
	std::vector< Eigen::Triplet< double > > triplets( numS * NumElementsPerSimplex * NumElementsPerSimplex );
	ThreadPool::ParallelFor
	(
		0 , numS ,
		[&]( size_t sIdx )
		{
			size_t idx = sIdx * NumElementsPerSimplex * NumElementsPerSimplex;
			for( unsigned int e1=0 ; e1<NumElementsPerSimplex ; e1++ ) for( unsigned int e2=0 ; e2<NumElementsPerSimplex ; e2++ , idx++ )
				triplets[ idx ] = Eigen::Triplet< double >( Idx(sIdx,e1) , Idx(sIdx,e2) , 0. );
		}
	);
	_M->setFromTriplets( triplets.begin() , triplets.end() );

	std::vector< std::map< size_t , double * > > entries( numF );
	// Using the fact that the matrix is column major
	ThreadPool::ParallelFor
	(
		0 , _M->outerSize() ,
		[&]( size_t c ){ for( typename Eigen::SparseMatrix< double >::InnerIterator it( *_M , c ) ; it ; ++it ) entries[c][ it.row() ] = &it.valueRef(); }
	);

	ThreadPool::ParallelFor
	(
		0 , numS ,
		[&]( size_t sIdx )
		{
			size_t idx = sIdx * NumElementsPerSimplex * NumElementsPerSimplex;
			for( unsigned int e1=0 ; e1<NumElementsPerSimplex ; e1++ ) for( unsigned int e2=0 ; e2<NumElementsPerSimplex ; e2++ , idx++ )
				_matrixEntries[ idx ] = entries[ Idx(sIdx,e1) ][ Idx(sIdx,e2) ];
		}
	);
}

template< unsigned int NumElementsPerSimplex >
void EigenMatrixEntries< NumElementsPerSimplex >::clear( void )
{
	ThreadPool::ParallelFor
	(
		0 , _M->outerSize() ,
		[&]( size_t o ){ for( typename Eigen::SparseMatrix< double >::InnerIterator it( *_M , o ) ; it ; ++it ) it.valueRef() = 0; }
	);
}