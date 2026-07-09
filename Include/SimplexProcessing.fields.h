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

		//////////////////////////////////////////////////////////////////////////////
		// A wrapper adding functionality to return the differential, as a function //
		//////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , typename T , typename Field >
		struct DifferentialFieldWrapper
		{
			struct Differential
			{
				Differential( const Field & f ) requires HasSimplexFunctionDifferential< Field , K , T > : _f(f){}
				SimplexProcessing::Differential< K , T > operator()( Position< K > p ) const { return _f.d(p); }
			protected:
				const Field _f;
			};

			// Return a function evaluating the differential
			Differential differential( void ) const requires HasSimplexFunctionDifferential< Field , K , T > { return Differential( static_cast< const Field & >( *this ) ); }
		};

		////////////////////////////////////////////////////
		// A structure for extracting the norm of a field //
		// [NOTE] The field is stored by value            //
		////////////////////////////////////////////////////
		template< unsigned int K , HasDotProduct T , HasSimplexFunction< K , T > Field >
		struct NormalizationField : public DifferentialFieldWrapper< K , T , NormalizationField< K , T , Field > >
		{
			static T Value( const Field & f , Position< K > p ){ return NormalizedValue( f(p) ); }
			static Differential< K , T > DValue( const Field & f , Position< K > p  ) requires HasSimplexFunctionDifferential< Field , K , T > { return DNormalizedValue( f(p) , f.d(p) ); }

			// Constructor
			NormalizationField( const Field & f ) : _f(f){}

			// Value evaluation
			T operator()( Position< K > p ) const { return Value( _f , p ); }

			// Derivative evaluation
			Differential< K , T > d( Position< K > p ) const requires HasSimplexFunctionDifferential< Field , K , T > { return DValue( _f , p ); }

		protected:
			Field _f;
		};

		// HELPER CLASSES 
		/////////////////

		//////////////////////////////
		// BASIC FUNCTIONS ON THE MESH

		///////////////////////////////////////////////////////////////
		// A structure for evaluating a linear function on a simplex //
		///////////////////////////////////////////////////////////////
		template< unsigned int K , typename _T >
		struct LinearInterpolant : public DifferentialFieldWrapper< K , _T , LinearInterpolant< K , _T > >
		{
			using T = _T;

			// [Class members]
			static T Value( const T x[K+1] , Position< K > p );
			static Differential< K , T > DValue( const T x[K+1] , Position< K > p );

			// The constructors
			LinearInterpolant( void );
			LinearInterpolant( const T x[K+1] );

			// The evaluation of the linear interpolation at a barycentric coordinate
			T operator()( Position< K > p ) const { return Value( _x , p ); }

			// The evaluation differential of the linear interpolant
			Differential< K , T > d( Position< K > p ) const { return DValue( _x , p ); }

			// Returns access to the coefficient at the prescribed vertex
			double &operator[]( unsigned int k ){ return _x[k]; }
			const double &operator[]( unsigned int k ) const { return _x[k]; }

		protected:
			T _x[K+1];
		};

		////////////////////////////
		// Rodrigues vector field //
		////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		struct PhongRodriguesVectorField : public DifferentialFieldWrapper< K , Point< double , N > , PhongRodriguesVectorField< K , N , Modulate > >
		{
			using T = Point< double , N >;

			// [Class members]

			static T Value( const T n[K+1] , const T x[K+1] , Position< K > p );
			static Differential< K , T > DValue( const T n[K+1] , const T x[K+1] , Position< K > p );

			// [Object members]

			PhongRodriguesVectorField( void ){}
			PhongRodriguesVectorField( const T n[K+1] ){ for( unsigned int k=0 ; k<=K ; k++ ) _n[k] = n[k]; }
			PhongRodriguesVectorField( const T n[K+1] , const T x[K+1] ) : PhongRodriguesVectorField( n ) { for( unsigned int k=0 ; k<=K ; k++ ) _x[k] = x[k]; }
			PhongRodriguesVectorField( const Simplex< double , N , K > & n , const T x[K+1] ) : PhongRodriguesVectorField( &n[0] , x ){}

			T &operator[]( unsigned int k ){ return _x[k]; }
			const T &operator[]( unsigned int k ) const { return _x[k]; }

			// The evaluation of the vector field
			T operator()( Position< K > p ) const { return Value( _n , _x , p ); }

			// The evaluation of the differential of the vector field
			Differential< K , T > d( Position< K > p ) const { return DValue( _n , _x , p ); }

		protected:
			static auto /* = std::pair< std::function< T (T) > , std::function< T (T) > > */ _DTransforms( Position< K > p , const T n[K+1] , const T x[K+1] );

			T _n[K+1] , _x[K+1];
		};

		// BASIC FUNCTIONS ON THE MESH
		//////////////////////////////

		/////////////////////////////////
		// FUNCTIONS OF THE MESH GEOMETRY

		///////////////////////////////////////////////////////////////
		// A field giving the (pseudo) differential of the embedding //
		///////////////////////////////////////////////////////////////
		template< unsigned int K >
		struct PhongRodriguesIntrinsicToExtrinsicTangentXFormField : public DifferentialFieldWrapper< K , Matrix< double , K , K+1 > , PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K > >
		{
			static const unsigned int Dim = K+1;
			using T = Matrix< double , K , Dim >;

			static T Value( const Point< double , Dim > & normal , const Point< double , Dim > normals[K+1] , const Matrix< double , K , Dim > & xForm , Position< K > p );
			static Differential< K , T > DValue( const Point< double , Dim > & normal , const Point< double , Dim > normals[K+1] , const Matrix< double , K , Dim > & xForm , Position< K > p );

			PhongRodriguesIntrinsicToExtrinsicTangentXFormField( void ){}
			PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			PhongRodriguesIntrinsicToExtrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );

			T operator()( Position< K > p ) const { return Value( _normal , _normals , _xForm , p ); }
			Differential< K , T > d( Position< K > p ) const { return DValue( _normal , _normals , _xForm , p ); }

		protected:
			Point< double , Dim > _normal , _normals[K+1];
			Matrix< double , K , Dim > _xForm;
		};

		//////////////////////////////////////////////////////////////////////////////
		// A field giving the inverse of the (pseudo) differential of the embedding //
		//////////////////////////////////////////////////////////////////////////////
		template< unsigned int K >
		struct PhongRodriguesExtrinsicToIntrinsicTangentXFormField
			: public DifferentialFieldWrapper< K , Matrix< double , K+1 , K > , PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > >
			, public PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >
		{
			static const unsigned int Dim = K+1;
			using T = Matrix< double , Dim , K >;

			using Differential = DifferentialFieldWrapper< K , Matrix< double , K+1 , K > , PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > >::Differential;
			using DifferentialFieldWrapper< K , Matrix< double , K+1 , K > , PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K > >::differential;

			static T Value( const SquareMatrix< double , K > & gInv , const Matrix< double , K , Dim > & i2e );
			static SimplexProcessing::Differential< K , T > DValue( const SquareMatrix< double , K > & gInv , const SimplexProcessing::Differential< K , Matrix< double , K , Dim > > & di2e );

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , Dim > & normal , const Point< double , Dim > normals[K+1] , const Matrix< double , K , Dim > & xForm , Position< K > p );
			static SimplexProcessing::Differential< K , T > DValue( const SquareMatrix< double , K > & gInv , const Point< double , Dim > & normal , const Point< double , Dim > normals[K+1] , const Matrix< double , K , Dim > & xForm , Position< K > p );

			PhongRodriguesExtrinsicToIntrinsicTangentXFormField( void ){}
			PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			PhongRodriguesExtrinsicToIntrinsicTangentXFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , p ); }
			SimplexProcessing::Differential< K , T > d( Position< K > p ) const{ return DValue( _gInv , _normal , _normals , _xForm , p ); }
		protected:
			SquareMatrix< double , K > _gInv;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_xForm;
		};

		//////////////////////////////////////////////////////////////////////////////////////////////////
		// A field giving the connection coefficients defined by (pseudo) differential of the embedding //
		//////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , bool Symmetrize=false >
		struct ConnectionCoefficientField : public PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >
		{
			static const unsigned int Dim = K+1;
			using T = AutoDiff::Tensor< K , K , K >;

			static T Value( const SquareMatrix< double , K > & gInv , const Matrix< double , K , Dim > & i2e , const Differential< K , Matrix< double , K , Dim > > & di2e );
			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , Dim > & normal , const Point< double , Dim > normals[K+1] , const Matrix< double , K , Dim > & xForm , Position< K > p );

			ConnectionCoefficientField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			ConnectionCoefficientField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , p ); }

		protected:
			SquareMatrix< double , K > _gInv;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_xForm;
		};

		///////////////////////////////////////
		// The second fundamental form field //
		///////////////////////////////////////
		template< unsigned int K , bool DifferentiateNormals=true >
		struct SecondFundamentalFormField : public PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >
		{
			static const unsigned int Dim = K+1;
			using T = SquareMatrix< double , K >;

			static T Value( const Point< double , Dim > & normal , const Point< double , Dim > normals[K+1] , const Matrix< double , K , Dim > & xForm , Position< K > p );

			SecondFundamentalFormField( const Point< double , Dim > vertices[K+1] , const Point< double , Dim > normals[K+1] );
			SecondFundamentalFormField( const Simplex< double , Dim , K > & vertices , const Simplex< double , Dim , K > & normals );

			// Given a simplex, returns a function returning the second fundamental form (expressed in the space of tangent vector fields) at any point in the simplex
			T operator()( Position< K > p ) const { return Value( _normal , _normals , _xForm , p ); }

		protected:
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesIntrinsicToExtrinsicTangentXFormField< K >::_xForm;
		};

		// FUNCTIONS OF THE MESH GEOMETRY
		/////////////////////////////////

		////////////////////////////////
		// GENERIC FUNCTIONS ON THE MESH

		////////////////////////////////////////////////////////////////////////////////////
		// Coefficients of a GENERIC vector-field vector field w.r.t. the coordinate axes //
		////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct IntrinsicVectorField
			: public DifferentialFieldWrapper< K , Point< double , K > , IntrinsicVectorField< K , N , VectorField > >
			, public PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >
		{
			using T = Point< double , K >;

			using Differential = DifferentialFieldWrapper< K , Point< double , K > , IntrinsicVectorField< K , N , VectorField > >::Differential;
			using DifferentialFieldWrapper< K , Point< double , K > , IntrinsicVectorField< K , N , VectorField > >::differential;

			static T Value( const Matrix< double , N , K > & e2i , const Point< double , N > & v );
			static SimplexProcessing::Differential< K , T > DValue( const Matrix< double , N , K > & e2i , const SimplexProcessing::Differential< K , Matrix< double , N , K > > & de2i , const Point< double , N > & v , const SimplexProcessing::Differential< K , Point< double , N > > & dv );

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const VectorField & vf , Position< K > p );
			static SimplexProcessing::Differential< K , T > DValue( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const VectorField & vf , Position< K > p );

			IntrinsicVectorField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & vf ) : PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >( vertices , normals ) , _vf( vf ) {}
			IntrinsicVectorField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & vf ) : IntrinsicVectorField( &vertices[0] , &normals[0] , vf ){}

			// The evaluation of the vector field
			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , _vf , p ); }

			// The evaluation of the differential of the vector field
			SimplexProcessing::Differential< K , T > d( Position< K > p ) const { return DValue( _gInv , _normal , _normals , _xForm , _vf , p ); }

		protected:
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_gInv;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_xForm;
			VectorField _vf;
		};

		////////////////////////////////////////////////////////////////////
		// The (intrinsic) covariant derivative of a GENERIC vector-field //
		////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct CovariantDerivativeField : public PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >
		{
			using T = SquareMatrix< double , K >;

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const VectorField & vf , Position< K > p );

			CovariantDerivativeField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & vf ) : PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >( vertices , normals ) , _vf( vf ) {}
			CovariantDerivativeField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & vf ) : CovariantDerivativeField( &vertices[0] , &normals[0] , vf ){}

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , _vf , p ); }

		protected:
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_gInv;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_xForm;
			VectorField _vf;
		};

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// The (extrinsic) covariant directional derivative of one GENERIC vector-field with respect to another //
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunction< K , Point< double , N > > DirectionField , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct CovariantDirectionalDerivativeField : public PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >
		{
			using T = Point< double , N >;

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const DirectionField & dir , const VectorField & vf , Position< K > p );

			CovariantDirectionalDerivativeField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const DirectionField & dir , const VectorField & vf ) : PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >( vertices , normals ) , _dir(dir) , _vf( vf ) {};
			CovariantDirectionalDerivativeField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const DirectionField & dir , const VectorField & vf ) : CovariantDirectionalDerivativeField( &vertices[0] , &normals[0] , dir , vf ){}

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , _dir , _vf , p ); }

		protected:
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_gInv;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_xForm;
			DirectionField _dir;
			VectorField _vf;
		};

		//////////////////////////////////////////////
		// The divergence of a GENERIC vector-field //
		//////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField >
		struct DivergenceField : public CovariantDerivativeField< K , N , VectorField >
		{
			using T = double;

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const VectorField & vf , Position< K > p );

			DivergenceField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField & vf ) : CovariantDerivativeField< K , N , VectorField >( vertices , normals , vf ) {};
			DivergenceField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField & vf ) : DivergenceField( &vertices[0] , &normals[0] , vf ){}

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , _vf , p ); }
		protected:
			using CovariantDerivativeField< K , N , VectorField >::_gInv;
			using CovariantDerivativeField< K , N , VectorField >::_normal;
			using CovariantDerivativeField< K , N , VectorField >::_normals;
			using CovariantDerivativeField< K , N , VectorField >::_xForm;
			using CovariantDerivativeField< K , N , VectorField >::_vf;
		};

		/////////////////////////////////////////////////////////////
		// The (extrinsic) commutator of two GENERIC vector-fields //
		/////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField1 , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField2 >
		struct CommutatorField : public PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >
		{
			using T = Point< double , N >;

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const VectorField1 & vf1 , const VectorField2 & vf2 , Position< K > p );

			CommutatorField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField1 & vf1 , const VectorField1 & vf2 ) : PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >( vertices , normals ) , _vf1(vf1) , _vf2(vf2) {};
			CommutatorField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField1 & vf1 , const VectorField1 & vf2 ) : CommutatorField( &vertices[0] , &normals[0] , vf1 , vf2 ){}

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , _vf1 , _vf2 , p ); }

		protected:
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_gInv;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_xForm;
			VectorField1 _vf1;
			VectorField2 _vf2;
		};

		//////////////////////////////////////////////////////////////
		// The (extrinsic) Lie bracket of two GENERIC vector-fields //
		//////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField1 , HasSimplexFunctionAndFunctionDifferential< K , Point< double , N > > VectorField2 >
		struct LieBracketField : public PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >
		{
			using T = Point< double , N >;

			static T Value( const SquareMatrix< double , K > & gInv , const Point< double , N > & normal , const Point< double , N > normals[K+1] , const Matrix< double , K , N > & xForm , const VectorField1 & vf1 , const VectorField2 & vf2 , Position< K > p );

			LieBracketField( const Point< double , N > vertices[K+1] , const Point< double , N > normals[K+1] , const VectorField1 & vf1 , const VectorField1 & vf2 ) : PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >( vertices , normals ) , _vf1(vf1) , _vf2(vf2) {};
			LieBracketField( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const VectorField1 & vf1 , const VectorField1 & vf2 ) : LieBracketField( &vertices[0] , &normals[0] , vf1 , vf2 ){}

			T operator()( Position< K > p ) const { return Value( _gInv , _normal , _normals , _xForm , _vf1 , _vf2 , p ); }

		protected:
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_gInv;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normal;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_normals;
			using PhongRodriguesExtrinsicToIntrinsicTangentXFormField< K >::_xForm;
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
		struct _SinglePhongRodriguesFunctionality : GenericVectorFieldFunctionality< K , N , PhongRodriguesVectorField< K , N , Modulate > >
		{
			using VectorField = PhongRodriguesVectorField< K , N , Modulate >;

			_SinglePhongRodriguesFunctionality( const Simplex< double , N , K > & vertices , const Simplex< double , N , K > & normals , const Point< double , N > vf[K+1] )
				requires std::constructible_from< GenericVectorFieldFunctionality< K , N , VectorField > , const Point< double , N > * , const Point< double , N > * , const VectorField & >
			: GenericVectorFieldFunctionality< K , N , VectorField >( &vertices[0] , &normals[0] , VectorField( &normals[0] , vf ) ) {}
		};

		template< template < unsigned int , unsigned int , typename , typename > typename GenericVectorFieldFunctionality , unsigned int K , unsigned int N , bool Modulate >
		struct _DoublePhongRodriguesFunctionality : GenericVectorFieldFunctionality< K , N , PhongRodriguesVectorField< K , N , Modulate > , PhongRodriguesVectorField< K , N , Modulate > >
		{
			using VectorField = PhongRodriguesVectorField< K , N , Modulate >;

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

		/////////////////////////////////////////////////////////////////////
		// The (extrinsic) commutator of two Phong-Rodrigues vector-fields //
		/////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesCommutatorField = _DoublePhongRodriguesFunctionality< CommutatorField , K , N , Modulate >;

		//////////////////////////////////////////////////////////////////////
		// The (extrinsic) Lie-Bracket of two Phong-Rodrigues vector-fields //
		//////////////////////////////////////////////////////////////////////
		template< unsigned int K , unsigned int N , bool Modulate=true >
		using PhongRodriguesLieBracketField = _DoublePhongRodriguesFunctionality< LieBracketField , K , N , Modulate >;

		// PHONG-RODRIGUES FUNCTIONS ON THE MESH
		////////////////////////////////////////

#include "SimplexProcessing.fields.inl"
	}
}
#endif // SIMPLEX_PROCESSING_FIELDS_INCLUDED
