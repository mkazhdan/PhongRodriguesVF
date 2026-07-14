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
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , HasElementIndexFunctor ElementIndex , HasMeshFunction< K , Point< double , NumElementsPerSimplex > > SystemVectorField >
Eigen::VectorXd RiemannianMesh< K , MeshType >::_systemVector( size_t fNum , ElementIndex && Idx , SystemVectorField && SVF , bool needsScaling ) const
{
	if( needsScaling )
	{
		auto _SVF = ScaledField< K >( measureScaleField() , std::forward< SystemVectorField >( SVF ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Vector< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , _SVF );
	}
	else
	{
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Vector< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , SVF );
	}
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , HasElementIndexFunctor ElementIndex , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > SystemField >
Eigen::SparseMatrix< double > RiemannianMesh< K , MeshType >::_systemMatrix( size_t fNum , ElementIndex && Idx , SystemField && Sys , bool needsScaling ) const
{
	if( needsScaling )
	{
		auto _Sys = ScaledField< K >( measureScaleField() , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , _Sys );
	}
	else
	{
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template Matrix< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) , Sys );
	}
}

template< unsigned int K , typename MeshType >
template< unsigned int NumElementsPerSimplex , HasElementIndexFunctor ElementIndex >
SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > RiemannianMesh< K , MeshType >::_eigenMatrixEntries( size_t fNum , ElementIndex && Idx ) const
{
	return SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex >( fNum , simplexNum() , std::forward< ElementIndex >( Idx ) );
}

template< unsigned int K , typename MeshType >
template< unsigned int QuadratureSamples , unsigned int NumElementsPerSimplex , HasMeshFunction< K , SquareMatrix< double , NumElementsPerSimplex > > SystemField >
void RiemannianMesh< K , MeshType >::_setSystemMatrixEntries( SystemIntegration::EigenMatrixEntries< NumElementsPerSimplex > & eme , SystemField && Sys , bool needsScaling ) const
{
	if( needsScaling )
	{
		auto _Sys = ScaledField< K >( measureScaleField() , std::forward< SystemField >( Sys ) );
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template SetMatrixEntries< NumElementsPerSimplex >( eme , simplexNum() , _Sys );
	}
	else
		return SystemIntegration::MCIntegrator< K , QuadratureSamples >::template SetMatrixEntries< NumElementsPerSimplex >( eme , simplexNum() , Sys );
}
