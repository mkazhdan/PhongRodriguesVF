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

////////////////////
// RiemannianMesh //
////////////////////
template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples >
double RiemannianMesh< K , MeshType >::measure( void ) const
{
	auto H = [&]( size_t sIdx ){ return [_S=measureScaleField(sIdx)]( SimplexProcessing::Position< K > p ){ return _S(p); }; };
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< double >( simplexNum() , SimplexProcessing::ArrayWrapper( H ) );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , typename T , HasMeshFunction< K , T > ScalarField >
T RiemannianMesh< K , MeshType >::integral( ScalarField && F ) const
{
	auto H = [&]( size_t sIdx ){ return [_F=F[sIdx],_S=measureScaleField(sIdx)]( SimplexProcessing::Position< K > p ){ return _F(p) * _S(p); }; };
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< T >( simplexNum() , SimplexProcessing::ArrayWrapper( H ) );
}

template< unsigned int K , typename MeshType >
template< typename T , HasMeshScaleFactorFunction< K > ScaleFactorField >
auto RiemannianMesh< K , MeshType >::_scaledIdentityField( ScaleFactorField && S ) const
{
	return SimplexProcessing::ArrayWrapper( [S]( size_t sIdx ){ return[s=S[sIdx]]( SimplexProcessing::Position< K > p ){ return [_s=s(p)]( T t ){ return t * _s; }; }; } );
}

template< unsigned int K , typename MeshType >
template< typename T , HasMeshScaleFactorFunction< K > ScaleFactorField >
auto RiemannianMesh< K , MeshType >::_scaledInverseMetricTensorField( ScaleFactorField && S ) const
{
	return SimplexProcessing::ArrayWrapper
	(
		[S,GInv=inverseMetricTensorField()]( size_t sIdx )
		{
			return [_GInv=GInv[sIdx],_S=S[sIdx]]( SimplexProcessing::Position< K > p )
			{
				return [s=_S(p),gInv=_GInv(p)]( SimplexProcessing::Differential< K , T > t ){ return gInv(t)*s; };
			};
		}
	);
}

template< unsigned int K , typename MeshType >
template< typename Elements >
auto RiemannianMesh< K , MeshType >::_DifferentialElements( Elements && E )
{
	return SimplexProcessing::ArrayWrapper
	(
		[&]( size_t sIdx )
		{
			return SimplexProcessing::ArrayWrapper
			(
				[_E=E[sIdx]]( size_t n )
				{
					return [__E=_E[n]]( SimplexProcessing::Position< K > p ){ return __E.d(p); };
				}
			);
		}
	);
}

template< unsigned int K , typename MeshType >
template< typename Elements >
auto RiemannianMesh< K , MeshType >::_ValueAndDifferentialElements( Elements && E )
{
	return SimplexProcessing::ArrayWrapper
	(
		[&]( size_t sIdx )
		{
			return SimplexProcessing::ArrayWrapper
			(
				[_E=E[sIdx]]( size_t n )
				{
					return [__E=_E[n]]( SimplexProcessing::Position< K > p )
					{
						return std::make_pair( __E(p) , __E.d(p) );
					};
				}
			);
		}
	);
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshFunction< K , T > ValueField , typename Elements , typename ElementIndex , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::VectorXd RiemannianMesh< K , MeshType >::_dual( size_t fNum , Elements && E , ValueField && F , ElementIndex && Idx , ScaleFactorField && S ) const
{
	auto VF = SystemIntegration::ScaledSystemVectorField< K , NumElementsPerSimplex , T >( std::forward< ScaleFactorField >( S ) , std::forward< ValueField >( F ) , std::forward< Elements >( E ) , this->template _scaledIdentityField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Vector< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , VF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshFunction< K , SimplexProcessing::Differential< K , T > > DifferentialField , typename Elements , typename ElementIndex , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::VectorXd RiemannianMesh< K , MeshType >::_dual( size_t fNum , Elements && E , DifferentialField && F , ElementIndex && Idx , ScaleFactorField && S ) const
{
	auto VF = SystemIntegration::ScaledSystemVectorField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( std::forward< ScaleFactorField >( S ) , std::forward< DifferentialField >( F ) , _DifferentialElements( E ) , this->template _scaledInverseMetricTensorField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Vector< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , VF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshFunction< K , T > ValueField , typename Elements , typename ElementIndex >
Eigen::VectorXd RiemannianMesh< K , MeshType >::_dual( size_t fNum , Elements && E , ValueField && F , ElementIndex && Idx ) const
{
	auto VF = SystemIntegration::SystemVectorField< K , NumElementsPerSimplex , T >( std::forward< ValueField >( F ) , std::forward< Elements >( E ) , this->template _scaledIdentityField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Vector< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , VF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , HasMeshFunction< K , SimplexProcessing::Differential< K , T > > DifferentialField , typename Elements , typename ElementIndex >
Eigen::VectorXd RiemannianMesh< K , MeshType >::_dual( size_t fNum , Elements && E , DifferentialField && F , ElementIndex && Idx ) const
{
	auto VF = SystemIntegration::SystemVectorField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( std::forward< DifferentialField >( F ) , _DifferentialElements( E ) , this->template _scaledInverseMetricTensorField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Vector< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , VF );
}


template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_mass( size_t fNum , Elements && E , ElementIndex && Idx ) const
{
	auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , T >( std::forward< Elements >( E ) , this->template _scaledIdentityField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_mass( size_t fNum , Elements && E , ElementIndex && Idx , ScaleFactorField && S ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , T >( std::forward< ScaleFactorField >( S ) , std::forward< Elements >( E ) , this->template _scaledIdentityField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_stiffness( size_t fNum , Elements && E , ElementIndex && Idx ) const
{
	auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( _DifferentialElements( E ) , this->template _scaledInverseMetricTensorField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex , HasMeshScaleFactorFunction< K > ScaleFactorField >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_stiffness( size_t fNum , Elements && E , ElementIndex && Idx , ScaleFactorField && S ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( std::forward< ScaleFactorField >( S ) , _DifferentialElements( E ) , this->template _scaledInverseMetricTensorField< T >( measureScaleField() ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
}


template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex , HasMeshSystemLinearMapOrBilinearFormFunction< K , NumElementsPerSimplex , T > SystemField >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_system( size_t fNum , Elements && E , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const
{
	if( needsScaling )
	{
		auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , T >( measureScaleField() , std::forward< Elements >( E ) , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
	}
	else
	{
		auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , T >( std::forward< Elements >( E ) , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
	}
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex , HasMeshDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_system( size_t fNum , Elements && E , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const
{
	if( needsScaling )
	{
		auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( measureScaleField() , _DifferentialElements( E ) , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
	}
	else
	{
		auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( _DifferentialElements( E ) , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , MF );
	}
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , typename ElementIndex , HasMeshSystemAndDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_system( size_t fNum , Elements && E , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const
{
	if( needsScaling )
	{
		auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , std::pair< T , SimplexProcessing::Differential< K , T > > >( measureScaleField() , _ValueAndDifferentialElements( E ) , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx )  , MF );
	}
	else
	{
		auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , std::pair< T , SimplexProcessing::Differential< K , T > > >( _ValueAndDifferentialElements( E ) , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx )  , MF );
	}
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexMass( size_t sIdx , Elements && E ) const
{
	auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , T >( std::forward< Elements >( E ) , this->template _scaledIdentityField< T >( measureScaleField() ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshScaleFactorFunction< K > ScaleFactorField >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexMass( size_t sIdx , Elements && E , ScaleFactorField && S ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , T >( std::forward< ScaleFactorField >( S ) , std::forward< Elements >( E ) , this->template _scaledIdentityField< T >( measureScaleField() ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexStiffness( size_t sIdx , Elements && E ) const
{
	auto MF = SystemIntegration::SystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( _DifferentialElements( E ) , this->template _scaledInverseMetricTensorField< T >( measureScaleField() ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshScaleFactorFunction< K > ScaleFactorField >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexStiffness( size_t sIdx , Elements && E , ScaleFactorField && S ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( std::forward< ScaleFactorField >( S ) , _DifferentialElements( E ) , this->template _scaledInverseMetricTensorField< T >( measureScaleField() ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshSystemLinearMapOrBilinearFormFunction< K , NumElementsPerSimplex , T > SystemField >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexSystem( size_t sIdx , Elements && E , SystemField && Sys ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , T >( measureScaleField() , std::forward< Elements >( E ) , std::forward< SystemField >( Sys ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexSystem( size_t sIdx , Elements && E , SystemField && Sys ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( measureScaleField() , _DifferentialElements( E ) , std::forward< SystemField >( Sys ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshSystemAndDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
SquareMatrix< double , NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_simplexSystem( size_t sIdx , Elements && E , SystemField && Sys ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , std::pair< T , SimplexProcessing::Differential< K , T > > >( measureScaleField() , _ValueAndDifferentialElements( E ) , std::forward< SystemField >( Sys ) );
	return SimplexProcessing::SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Integral< SquareMatrix< double , NumElementsPerSimplex > >( MF[sIdx] );
}

template< unsigned int K , typename MeshType >
template< unsigned int NumElementsPerSimplex , typename ElementIndex >
SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_eigenMatrixEntries( size_t fNum , ElementIndex && Idx ) const
{
	return SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshSystemLinearMapOrBilinearFormFunction< K , NumElementsPerSimplex , T > SystemField >
void RiemannianMesh< K , MeshType >::_setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , Elements && E , SystemField && Sys ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , T >( measureScaleField() , std::forward< Elements >( E ) , std::forward< SystemField >( Sys ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template SetMatrixEntries< NumElementsPerSimplex >( eme , simplexNum() , MF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
void RiemannianMesh< K , MeshType >::_setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , Elements && E , SystemField && Sys ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , SimplexProcessing::Differential< K , T > >( measureScaleField() , _DifferentialElements( E ) , std::forward< SystemField >( Sys ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template SetMatrixEntries< NumElementsPerSimplex >( eme , simplexNum() , MF );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , typename T , typename Elements , HasMeshSystemAndDSystemFunction< K , NumElementsPerSimplex , T > SystemField >
void RiemannianMesh< K , MeshType >::_setSystemEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , Elements && E , SystemField && Sys ) const
{
	auto MF = SystemIntegration::ScaledSystemMatrixField< K , NumElementsPerSimplex , std::pair< T , SimplexProcessing::Differential< K , T > > >( measureScaleField() , _ValueAndDifferentialElements( E ) , std::forward< SystemField >( Sys ) );
	return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template SetMatrixEntries< NumElementsPerSimplex >( eme , simplexNum() , MF );
}

template< unsigned int K , typename MeshType >
template< unsigned int NumElementsPerSimplex , typename T , HasMeshTangentVectorFunction< K > TangentVectorField >
auto RiemannianMesh< K , MeshType >::_DerivationSystemField( TangentVectorField && VF )
{
	auto Sys = [&]( size_t sIdx )
	{
		return [sVF=VF[sIdx]]( SimplexProcessing::Position< K > p )
		{
			return [vf=sVF(p)]( const std::pair< T , SimplexProcessing::Differential< K , T > > * valuesAndDValues  )
			{
				SquareMatrix< double , NumElementsPerSimplex > S;
				for( unsigned int m=0 ; m<NumElementsPerSimplex ; m++ )
				{
					T value{};
					for( unsigned int k=0 ; k<K ; k++ ) value += valuesAndDValues[m].second[k] * vf[k];
					for( unsigned int n=0 ; n<NumElementsPerSimplex ; n++ ) S(n,m) = SimplexProcessing::DotProduct( value , valuesAndDValues[n].first );
				}
				return S;
			};
		};
	};
	return SimplexProcessing::ArrayWrapper( Sys );
}
