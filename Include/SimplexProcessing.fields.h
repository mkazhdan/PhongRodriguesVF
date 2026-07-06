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

#ifndef SIMPLEX_PROCESSING_FIELDS_INCLUDED
#define SIMPLEX_PROCESSING_FIELDS_INCLUDED

namespace MishaK
{
	namespace SimplexProcessing
	{
		/////////////////
		// HELPER CLASSES 

		//////////////////////////////////////////////////////////////////////////////////////
		// A structure for extracting the differential of a field as a stand-alone function //
		// [NOTE] The field is stored by value                                              //
		//////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , typename T , HasSimplexFunctionDifferential< K , T > Field >
		struct DifferentialField
		{
			DifferentialField( const Field & f ) : _f(f){}
			Differential< K , T > operator()( Position< K > p ) const { return _f.d(p); }
		protected:
			Field _f;
		};

		//////////////////////////////////////////////////////////////////////////////
		// A wrapper adding functionality to return the differential, as a function //
		//////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , typename T , HasSimplexFunctionDifferential< K , T > Field >
		struct DifferentialFieldWrapper : public Field
		{
			// Inherit base constructors
			using Field::Field;

			// Ensure that evaluation is still supported
			using Field::d;

			// The differential function type
			using Differential = DifferentialField< K , T , Field >;

			// Return a function evaluating the differential
			DifferentialField< K , T , Field > d( void ) const { return Differential( static_cast< const Field & >( *this ) ); }
		};

		/////////////////////////////////////////////////////////////////////////////////////////////////////
		// Functionality for validating that analytic derivatives are approximated by discrete derivatives //
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , HasDotProduct T >
		struct DerivativeTester
		{
			template< HasSimplexFunctionAndFunctionDifferential< K , T > Field >
			static double SquareError( Position< K > p , const Field & field , double delta );

			template< HasSimplexFunctionAndFunctionDifferential< K , T > Field >
			static double SquareError( const Field & field , unsigned int testCount , double delta );
		};

		////////////////////////////////////////////////////
		// A structure for extracting the norm of a field //
		// [NOTE] The field is stored by value            //
		////////////////////////////////////////////////////
		template< unsigned int K , HasDotProduct T , HasSimplexFunction< K , T > Field >
		struct _NormalizationField
		{
			// Constructor
			_NormalizationField( const Field & f ) : _f(f){}

			// Value evaluation
			T operator()( Position< K > p ) const { T t = _f(p) ; return t / sqrt( DotProduct( t , t ) ); }

			// Derivative evaluation
			Differential< K , T > d( Position< K > p ) const;
		protected:
			Field _f;
		};

		template< unsigned int K , typename T , typename Field >
		using NormalizationField = DifferentialFieldWrapper< K , T , _NormalizationField< K , T , Field > >;

		// HELPER CLASSES 
		/////////////////

		//////////////////////////////
		// BASIC FUNCTIONS ON THE MESH

		///////////////////////////////////////////////////////////////
		// A structure for evaluating a linear function on a simplex //
		///////////////////////////////////////////////////////////////
		template< unsigned int K , typename _T >
		struct _LinearInterpolant
		{
			using T = _T;

			// [Class members]

			// The constructors
			_LinearInterpolant( void );
			_LinearInterpolant( const T x[K+1] );

			// The evaluation of the linear interpolation at a barycentric coordinate
			T operator()( Position< K > p ) const { return Value( p , _x ); }

			// The evaluation differential of the linear interpolant
			Differential< K , T > d( Position< K > p ) const { return DValue( p , _x ); }

			// Returns access to the coefficient at the prescribed vertex
			double &operator[]( unsigned int k ){ return _x[k]; }
			const double &operator[]( unsigned int k ) const { return _x[k]; }

			// [Class members]
			static T Value( Position< K > p , const T x[K+1] );
			static Differential< K , T > DValue( Position< K > p , const T x[K+1] );

		protected:
			T _x[K+1];
		};

		template< unsigned int K , typename T >
		using LinearInterpolant = DifferentialFieldWrapper< K , T , _LinearInterpolant< K , T > >;

		////////////////////////////
		// Rodrigues vector field //
		////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		struct _PhongRodriguesVectorField
		{
			using T = Point< double , N >;

			// [Class members]

			static T Value( Position< K > p , const T n[K+1] , const T x[K+1] );
			static Differential< K , T > DValue( Position< K > p , const T n[K+1] , const T x[K+1] );

			// [Object members]

			_PhongRodriguesVectorField( void ){}
			_PhongRodriguesVectorField( const T n[K+1] ){ for( unsigned int k=0 ; k<=K ; k++ ) _n[k] = n[k]; }
			_PhongRodriguesVectorField( const T n[K+1] , const T x[K+1] ) : _PhongRodriguesVectorField( n ) { for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = x[k]; }
			_PhongRodriguesVectorField( const Simplex< double , N , K > & n , const T x[K+1] ) : _PhongRodriguesVectorField( &n[0] , x ){}

			T &operator[]( unsigned int k ){ return _x[k]; }
			const T &operator[]( unsigned int k ) const { return _x[k]; }

			// The evaluation of the vector field
			T operator()( Position< K > p ) const { return Value( p , _n , _x ); }

			// The evaluation of the differential of the vector field
			Differential< K , T > d( Position< K > p ) const { return DValue( p , _n , _x ); }

		protected:
			static auto /* = std::pair< std::function< T (T) > , std::function< T (T) > > */ _DTransforms( Position< K > p , const T n[K+1] , const T x[K+1] );

			T _n[K+1] , _x[K+1];
		};

		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesVectorField = DifferentialFieldWrapper< K , typename _PhongRodriguesVectorField< K , N , Modulate >::T , _PhongRodriguesVectorField< K , N , Modulate > >;

		// BASIC FUNCTIONS ON THE MESH
		//////////////////////////////

		/////////////////////////////////
		// FUNCTIONS OF THE MESH GEOMETRY

		///////////////////////////////////////////////////////////////
		// A field giving the (pseudo) differential of the embedding //
		///////////////////////////////////////////////////////////////
		template< unsigned int K >
		struct _PhongRodriguesIntrinsicToExtrinsicTangentXFormField
		{
			static const unsigned int Dim = K+1;
			using T = Matrix< double , K , Dim >;

			_PhongRodriguesIntrinsicToExtrinsicTangentXFormField( void ){}
			_PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			_PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );
			T operator()( Position< K > p ) const;
			Differential< K , T > d( Position< K > p ) const;
		protected:
			LinearInterpolant< K , Point< double , Dim > > _N;
			Point< double , Dim > _normal;
			Matrix< double , K , Dim > _xForm;
		};
		template< unsigned int K >
		using PhongRodriguesIntrinsicToExtrinsicTangentXFormField = DifferentialFieldWrapper< K , typename _PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::T , _PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K > >;

		//////////////////////////////////////////////////////////////////////////////
		// A field giving the inverse of the (pseudo) differential of the embedding //
		//////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		struct _PhongRodriguesExtrinsicToIntrinsicTangentXFormField
		{
			static const unsigned int Dim = K+1;
			using T = Matrix< double , Dim , K >;

			_PhongRodriguesExtrinsicToIntrinsicTangentXFormField( void ){}
			_PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			_PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );
			T operator()( Position< K > p ) const;
			Differential< K , T > d( Position< K > p ) const;
		protected:
			_PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >  _i2e;
			SquareMatrix< double , K > _gInv;
		};
		template< unsigned int K >
		using PhongRodriguesExtrinsicToIntrinsicTangentXFormField = DifferentialFieldWrapper< K , typename _PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::T , _PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > >;

		//////////////////////////////////////////////////////////////////////////////////////////////////
		// A field giving the connection coefficients defined by (pseudo) differential of the embedding //
		//////////////////////////////////////////////////////////////////////////////////////////////////

		template< unsigned int K >
		struct ConnectionCoefficientField
		{
			static const unsigned int Dim = K+1;
			using T = AutoDiff::Tensor< K , K , K >;

			ConnectionCoefficientField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			ConnectionCoefficientField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );
			T operator()( Position< K > p ) const;

		protected:
			PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K > _i2e;
			SquareMatrix< double , K > _gInv;
		};

		///////////////////////////////////////
		// The second fundamental form field //
		///////////////////////////////////////
		template< unsigned int K , bool DifferentiateNormals=true >
		struct SecondFundamentalFormField
		{
			static const unsigned int Dim = K+1;
			using T = SquareMatrix< double , K >;

			SecondFundamentalFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			SecondFundamentalFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );

			// Given a simplex, returns a function returning the second fundamental form (expressed in the space of tangent vector fields) at any point in the simplex
			T operator()( Position< K > p ) const;

		protected:
			PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K > _i2e;
			NormalizationField< K , Point< double , Dim > , LinearInterpolant< K , Point< double , Dim > > > _normals;
		};

		// FUNCTIONS OF THE MESH GEOMETRY
		/////////////////////////////////

		////////////////////////////////
		// GENERIC FUNCTIONS ON THE MESH

		////////////////////////////////////////////////////////////////////////////////////
		// Coefficients of a GENERIC vector-field vector field w.r.t. the coordinate axes //
		////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct _IntrinsicVectorField
		{
			using T = Point< double , K >;

			_IntrinsicVectorField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & vf ) : _e2i( vertices , normals ) , _vf( vf ) {}
			_IntrinsicVectorField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & vf ) : _IntrinsicVectorField( &vertices[0] , &normals[0] , vf ){}

			// The evaluation of the vector field
			T operator()( Position< K > p ) const;

			// The evaluation of the differential of the vector field
			Differential< K , T > d( Position< K > p ) const;

		protected:
			PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > _e2i;
			VectorField _vf;
		};

		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		using IntrinsicVectorField = DifferentialFieldWrapper< K , typename _IntrinsicVectorField< K , N , VectorField >::T , _IntrinsicVectorField< K , N , VectorField > >;

		////////////////////////////////////////////////////////////////////
		// The (intrinsic) covariant derivative of a GENERIC vector-field //
		////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct CovariantDerivativeField
		{
			using T = SquareMatrix< double , K >;

			CovariantDerivativeField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & vf ) : _e2i( vertices , normals ) , _vf( vf ) {}
			CovariantDerivativeField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & vf ) : CovariantDerivativeField( &vertices[0] , &normals[0] , vf ){}
			T operator()( Position< K > p ) const;
		protected:
			PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > _e2i;
			VectorField _vf;
		};

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// The (extrinsic) covariant directional derivative of one GENERIC vector-field with respect to another //
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunction< K , Point< double , N > > DirectionField , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct CovariantDirectionalDerivativeField
		{
			using T = Point< double , N >;

			CovariantDirectionalDerivativeField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const DirectionField & dir , const VectorField & vf ) : _e2i( vertices , normals ) , _dir(dir) , _vf( vf ) {};
			CovariantDirectionalDerivativeField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const DirectionField & dir , const VectorField & vf ) : CovariantDirectionalDerivativeField( &vertices[0] , &normals[0] , dir , vf ){}
			T operator()( Position< K > p ) const;
		protected:
			PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > _e2i;
			DirectionField _dir;
			VectorField _vf;
		};

		//////////////////////////////////////////////
		// The divergence of a GENERIC vector-field //
		//////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct DivergenceField
		{
			using T = double;

			DivergenceField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & vf ) : _dvf( vertices , normals , vf ) {};
			DivergenceField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & vf ) : DivergenceField( &vertices[0] , &normals[0] , vf ){}
			T operator()( Position< K > p ) const;
		protected:
			CovariantDerivativeField< K , N , VectorField > _dvf;
		};

		//////////////////////////////////////////////////////////
		// The (extrinsic) bracket of two GENERIC vector-fields //
		//////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField1 , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField2 >
		struct BracketField
		{
			using T = Point< double , N >;

			BracketField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField1 & vf1 , const VectorField1 & vf2 ) : _e2i( vertices , normals ) , _vf1(vf1) , _vf2(vf2) {};
			BracketField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField1 & vf1 , const VectorField1 & vf2 ) : BracketField( &vertices[0] , &normals[0] , vf1 , vf2 ){}
			T operator()( Position< K > p ) const;
		protected:
			PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > _e2i;
			VectorField1 _vf1;
			VectorField2 _vf2;
		};

		// GENERIC FUNCTIONS ON THE MESH
		////////////////////////////////


		////////////////////////////////////////
		// PHONG-RODRIGUES FUNCTIONS ON THE MESH
		
		/////////////////////////////////////////////////////////////////////////////////////
		// Wrappers for implementing GENERIC functionality using the Phong-Rodrigues basis //
		/////////////////////////////////////////////////////////////////////////////////////

		template< template < unsigned int , unsigned int , typename > typename GenericVectorFieldFunctionality , unsigned int K , unsigned int N , bool Modulate >
		struct _SinglePhongRodriguesFunctionality : GenericVectorFieldFunctionality< K , N , _PhongRodriguesVectorField< K , N , Modulate > >
		{
			using VectorField = _PhongRodriguesVectorField< K , N , Modulate >;

			_SinglePhongRodriguesFunctionality( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const Point< double , N > vf[K+1] )
				requires std::constructible_from< GenericVectorFieldFunctionality< K , N , VectorField > , const Point< double , N > * , const Point< double , N > * , const VectorField & >
			: GenericVectorFieldFunctionality< K , N , VectorField >( &vertices[0] , &normals[0] , VectorField( &normals[0] , vf ) ) {}
		};

		template< template < unsigned int , unsigned int , typename , typename > typename GenericVectorFieldFunctionality , unsigned int K , unsigned int N , bool Modulate >
		struct _DoublePhongRodriguesFunctionality : GenericVectorFieldFunctionality< K , N , _PhongRodriguesVectorField< K , N , Modulate > , _PhongRodriguesVectorField< K , N , Modulate > >
		{
			using VectorField = _PhongRodriguesVectorField< K , N , Modulate >;

			_DoublePhongRodriguesFunctionality( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const Point< double , N > vf1[K+1] , const Point< double , N > vf2[K+1] )
				requires std::constructible_from< GenericVectorFieldFunctionality< K , N , VectorField , VectorField > , const Point< double , N > * , const Point< double , N > * , const VectorField & , const VectorField & >
			: GenericVectorFieldFunctionality< K , N , VectorField , VectorField >( &vertices[0] , &normals[0] , VectorField( &normals[0] , vf1 ) , VectorField( &normals[0] , vf2 ) ) {}
		};

		///////////////////////////////////////////////////////////////////////////////
		// Coefficients of a Phong-Rodrigues vector-field w.r.t. the coordinate axes //
		///////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesIntrinsicVectorField = _SinglePhongRodriguesFunctionality< IntrinsicVectorField , K , N , Modulate >;

		////////////////////////////////////////////////////////////////////////////
		// The (intrinsic) covariant derivative of a Phong-Rodrigues vector-field //
		////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesCovariantDerivativeField = _SinglePhongRodriguesFunctionality< CovariantDerivativeField , K , N , Modulate >;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// The (extrinsic) covariant directional derivative of one Phong-Rodrigues vector-field with respect to another //
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesCovariantDirectionalDerivativeField = _DoublePhongRodriguesFunctionality< CovariantDirectionalDerivativeField , K , N , Modulate >;

		////////////////////////////////////////////////////////
		// The divergence of the Phong-Rodrigues vector-field //
		////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesDivergenceField = _SinglePhongRodriguesFunctionality< DivergenceField , K , N , Modulate >;

		//////////////////////////////////////////////////////////////////
		// The (extrinsic) bracket of two Phong-Rodrigues vector-fields //
		//////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesBracketField = _DoublePhongRodriguesFunctionality< BracketField , K , N , Modulate >;

		// PHONG-RODRIGUES FUNCTIONS ON THE MESH
		////////////////////////////////////////

#include "SimplexProcessing.fields.inl"
	}
}
#endif // SIMPLEX_PROCESSING_FIELDS_INCLUDED
