<center><h2>Phong-Rodrigues Extrinsic Vector-Field Processing (Version 1.15)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<a href="#EXAMPLES">examples</a>
<a href="#COMPILATION">compilation</a>
<a href="#USAGE">usage</a>
<a href="#CHANGES">changes</a>
</center>
<hr>
This software supports vector-field processing using the extrinsically defined Phong-Rodrigues basis. Supported applications include computation of:
<UL>
<LI>sparse interpolation
<LI>bracket of two vector-fields
</UL>
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
<b>Papers:</b>
<a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP26.pdf">[Liu, Stein, Vaxman, Ben-Chen, and Kazhdan, 2026]</a>
<br>
<b>Executables: </b>
<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.15/PRVF.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.15/PRVF.Source.zip">ZIP</a> <a href="https://github.com/mkazhdan/PhongRodriguesVF">GitHub</a><br>
<B>Data:</B>
<A HREF="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/PRVF.Data.zip">ZIP</A><br>
<b>Older Versions:</b>
<a href="https://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version1.10/">V1.10</a>,<BR>
<a href="https://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version1.00/">V1.00</a>
</ul>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>
<ul>

<dl>
<details>
<summary>
<font size="+1"><b>SparseInterpolation</b></font>:
Supports the construction of a tangent vector-field that interpolates the prescribed constraints while being as-smooth-as-possible everywhere else.<BR>
"Smoothness" is defined in terms of the connection, the Hodge, or the Killing energy.
</summary>
<dt><b>--in</b> &lt;<i>input mesh and interpolation constraints</i>&gt;</dt>
<dd> These two strings specify the the names of the files describing the mesh and the interpolation constraints.<BR>
The mesh is assumed to be in <A href="https://sites.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.<BR>
The constraints are assumed to be in ASCII format, with each constraint described by a white-space separated quadruple of rnumbers. The first is an integer value describing the index of the constrained vertex. The remaining three are floating points giving the extrinsic coordinates of the tangent vector-field at that vertex.<BR>
<B>Note</B>:
<UL>
<LI> If they are not provided, the executable will generate per-vertex normals.
<LI> Normals are rescaled to have unit-lengths.
<LI> Input tangent vector-fields are processed to have the normal components projected out.
</UL>
</dd>

<dt>[<b>--out</b> &lt;<i>output mesh</i>&gt;]</dt>
<dd> This string is the name of the file to which the interpolating vector-field will be written.</B>
If the file-name ends with the <I>.ply</I> extension, the output mesh and vector-field will be written in <A href="https://sites.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format with the vertex positions encoded into parmaeters <I>x</I>, <I>y</I>, and <I>z</I>, the normals encoded into parameters <I>nx</I>, <I>ny</I>, and <I>nz</I>, and the vector-field values encoded into parameters <I>vf0</I>, <I>vf1</I>, <I>vf2</I>.<BR>
Otherwise, the output vector-field will be written in ASCII format, with each vertex's tangent vector described by a white-space separated triple of floating-point values giving the extrinsic coordinates of the tangent vector-field at that vertex.
</dd>

<dt>[<b>--qSamples</b> &lt;<i>quadrature samples</i>&gt;]</dt>
<dd> This integer specifies the number of quadrature points used for estimating an integral over a triangle.<BR>
Supported values are in {1, 3, 4, 6, 7, 12, 13, 24, 27, 32}.<BR>
 The default value for this parameter is 3.
</dd>

<dt>[<b>--energy</b> &lt;<i>type of vector-field energy</i>&gt;]</dt>
<dd> This integer specifies the type of energy used to define "smoothness". Valid values are:
<UL>
<LI><b>0</B>: Connection
<LI><b>1</B>: Hodge
<LI><b>2</B>: Killing
</UL>
The default value for this parameter is 0.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If this flag is enabled, performance information is printed to <CODE>stdout</CODE>.

</details>
</dl>


<dl>
<details>
<summary>
<font size="+1"><b>Bracket</b></font>:
Supports the evaluation and computation of the bracked of two vector-fields.
</summary>
<dt><b>--in</b> &lt;<i>input mesh, vector-field 1, and vector-field 2</i>&gt;</dt>
<dd> These three strings specify the the names of the files describing the mesh and the two vector-fields whose bracket is to be computed.<BR>
The mesh is assumed to be in <A href="https://sites.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.<BR>
If the vector-field file-name ends with the <I>.ply</I> extension, the vector-field is assumed to be written in <A href="https://sites.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format with the vector-field values encoded into parameters <I>vf0</I>, <I>vf1</I>, <I>vf2</I>.<BR>
Otherwise, the vector-field values are assumed to be in ASCII format, written in the order of the vertices, with each vertex's vector-field value expressed a white-space separated triple of floating point values giving the extrinsic coordinates of the tangent vector-field at that vertex.<BR>
<B>Note</B>:
<UL>
<LI> If they are not provided, the executable will generate per-vertex normals.
<LI> Normals are rescaled to have unit-lengths.
<LI> Input tangent vector-fields are processed to have the normal components projected out.
</UL>
</dd>

<dt>[<b>--out</b> &lt;<i>output mesh</i>&gt;]</dt>
<dd> This string is the name of the file to which the bracket of the two vector-fields will be written.</B>
If the file-name ends with the <I>.ply</I> extension, the output mesh and vector-field will be written in <A href="https://sites.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format with the vertex positions encoded into parmaeters <I>x</I>, <I>y</I>, and <I>z</I>, the normals encoded into parameters <I>nx</I>, <I>ny</I>, and <I>nz</I>, and the vector-field values encoded into parameters <I>vf0</I>, <I>vf1</I>, <I>vf2</I>.<BR>
Otherwise, the output vector-field will be written in ASCII format, with each vertex's tangent vector described by a white-space separated triple of floating-point values giving the extrinsic coordinates of the tangent vector-field at that vertex.
</dd>

<dt>[<b>--qSamples</b> &lt;<i>quadrature samples</i>&gt;]</dt>
<dd> This integer specifies the number of quadrature points used for estimating an integral over a triangle.<BR>
Supported values are in {1, 3, 4, 6, 7, 12, 13, 24, 27, 32}.<BR>
 The default value for this parameter is 3.
</dd>

<dt>[<b>--covDiff</b>]</dt>
<dd> If this flag is enabled, the method computes the covariant derivatives of the two vector-fields, evaluates each of those along the other vector-field, and returns the difference.

<dt>[<b>--verbose</b>]</dt>
<dd> If this flag is enabled, performance information is printed to <CODE>stdout</CODE>.

</details>
</dl>


</ul>


<hr>
<a name="EXAMPLES"><b>EXAMPLES (WITH SAMPLE DATA)</b></a><br>
For testing purposes, <A HREF="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/PRVF.Data.zip">the provided archive</A> contains the "kitten" model with two tangent vector-field constraints, as well as two tangent vector-fields defined over the torus.

<ul>

<dl>
<details>
<summary>
<font size="+1"><b>SparseInterpolation</b></font>
</summary>
To run this executable you must specify the input mesh as well as the constraints:
<blockquote><code>% Bin/*/SparseInterpolation --in PRVF.Data/kitten.ply PRVF.Data/kitten.txt --out connection.ply</code></blockquote>
This computes the tangent vector-field with minimum connection energy that interpolates the two vertex constraints, and writes out the vector-field to <CODE>connection.ply</CODE>.<BR>
You can minimize the Hodge energy instead by running:
<blockquote><code>% Bin/*/SparseInterpolation --in PRVF.Data/kitten.ply PRVF.Data/kitten.txt --energy --out hodge.ply</code></blockquote>
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>Bracket</b></font>
</summary>
To run this executable you must specify the input mesh as well as the two vector-fields:
<blockquote><code>% Bin/*/SparseInterpolation --in PRVF.Data/torus.X.ply PRVF.Data/torus.X.ply torus.Y.ply --out torus.XY.ply</code></blockquote>
This evaluates the bracket and computes its best-fit representation in terms of the Phong-Rodriguess vector-field basis, writing out the mesh and vector-field to the file <CODE>torus.XY.ply</CODE>.<BR>
(Note that as the files <CODE>PRVF.Data/torus.X.ply</CODE> and <CODE>PRVF.Data/torus.Y.ply</CODE> encode both the geometry and the vector-fields, we pass in <CODE>PRVF.Data/torus.X.ply</CODE> twice, once for the endoding of the geometry and once for the vector-field.)
</details>
</dl>



</ul>


<hr>
<details>
<summary>
<a name="COMPILATION"><b>COMPILATION</b></a><br>
</summary>
<UL>
<LI>Compilation requires a sparse matrix library, as well as a linear solver. By default, we use the <CODE>LDLt</CODE> implementation provided by <A HREF="https://eigen.tuxfamily.org/">Eigen</A>. </UL>
</details>

<hr>
<details>
<summary>
<a name="USAGE"><b>USAGE</b></a><br>
</summary>
<OL>
<DL>

<DT><B>Set-up</B>
<DD>
The implementation is header-only.
<UL>
<LI> To perform Phong-Rodrigues vector-field processing, you will need to include the file <code>Include/SimplicialMeshProcessing.h</code>, which defines the functionality within the <code>MishaK::SimplicialMesh</code> namespace.
<LI>
Most of the functionality is defined by the class <CODE>EmbeddedPhongMesh&lt;2&gt;</CODE>.
<LI> To define a mesh object, call the constructor:
<BLOCKQUOTE><CODE>
EmbeddedPhongMesh&lt;2&gt;::EmbeddedPhongMesh( std::vector&lt; MishaK::Point&lt;double,3&gt; &gt; &amp; vertices , std::vector&lt; MishaK::Point&lt;double,3&gt; &gt; &amp; normals , std::vector&lt; MishaK::SimplexIndex&lt;2&gt; &gt; &amp; triangles );
</CODE></BLOCKQUOTE>
The geometry is represented as an indexed mesh, with:
<UL>
<LI>
The (shared) vertices and normals described using <code>std::vector</code>s of type <code>MishaK::Point&lt;double,3&gt;</code>, which acts as an array of size <code>3</CODE> and additionally supports algebraic operations.
<LI>The triangle connectivity is described using an <code>std::vector</code> of type <code>MishaK::SimplexIndex&lt;2&gt;</code>, which acts as an array of size <code>3</CODE>, giving the indices of the three corners of the triangle.
</UL>
Note that the constructed mesh captures the vertex, normal, and triangle information <b>by reference</b>, so geometry should not be deleted until after all the processing is completed.
</UL>

<DT><B>System Representation</B></DT>
<DD>
<UL>
<LI>System matrices are represented using <A HREF="https://libeigen.gitlab.io/eigen/docs-5.0/">Eigen</A>'s sparse-matrix <code>Eigen::SparseMatrix&lt;double&gt;</code> class and system vectors are represented using the <code>Eigen::VectorXd</code> class. For a triangle mesh with <code>|V|</CODE> vertices, the extrinsic representation of the vector-field has <code>3&middot;|V|</code> degrees of freedom, with the three Euclidean coordinates of the vector at vertex <code>i</code> encoded in entries <code>3&middot;i</code>, <code>3&middot;i+1</code>, and <code>3&middot;i+2</code>.
<LI> In practice, processing is done over <b>tangent</b> vector-fields, so that there are only two degrees of freedom per vertex. This can be realized by selecting a frame for each vertex, perpendicular to the vertex's normal, and describing the vertex's tangent vector with respect to the frame. To facilitate this, the class defines the member function:
<BLOCKQUOTE><CODE>
Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::tangentProlongation( void )
</CODE></BLOCKQUOTE>
which returns a sparse <code>3&middot;|V|&times;2&middot;|V|</code> matrix mapping a representation of tangent vectors in terms of two degrees of freedom at a vertex (with respect to a frame perpendicular to the normal) to the Euclidean coordinates of the vectors at the vertices. The frame is chosen to be orthonormal, so that if <CODE>P</CODE> is the matrix returned by the method, the matrix <CODE>P.transpose()*P</CODE> is the <code>3&middot;|V|&times;3&middot;|V|</code> matrix describing the operation of projecting out the normal component. Thus, for example:
<UL>
<LI>Given a symmetric <code>3&middot;|V|&times;2&middot;|V|</code> matrix <CODE>Q</CODE> describing the energy of a vector-field expressed in terms Euclidean coordinates, the matrix representing the restriction of the energy to tangent-vectors is:
<BLOCKQUOTE><CODE>P.transpose()*Q*P</CODE></BLOCKQUOTE>
<LI>Given a <code>2&middot;|V|</code>-dimensional vector <code>v</code> describing the tangent-vectors at the vertices, with respect to the frame, the representation of the tangent-vectors in terms of Euclidean coordinates is:
<BLOCKQUOTE><CODE>P*v</CODE></BLOCKQUOTE>
<LI>Given a <code>3&middot;|V|</code>-dimensional vector <code>d</code> describing the dual of a vector-field with respect to the Euclidean Phong-Rodrigues basis (i.e. the integral of the vector-field against each of the <code>3&middot;|V|</code> basis vector-fields), the dual representation with respect to the tangent Phong-Rodrigues basis is:
<BLOCKQUOTE><CODE>P.transpose()*d</CODE></BLOCKQUOTE>
<LI>Finally, because the implementation defines the prolongation matrix <code>P</code> using an orthonormal frame, given a <code>3&middot;|V|</code>-dimensional vector <code>v</code> describing the tangent-vectors at the vertices in Euclidean coordinates the representation of the tangent-vectors with respect to the frame is:
<BLOCKQUOTE><CODE>P*v</CODE></BLOCKQUOTE>
</UL>
</UL>

<DT><B>System Energies</B>
<DD> The system supports computing symmetric, positive (semi-)definite <code>3&middot;|V|&times;3&middot;|V|</code> matrices representing different energies. All the functionalities are templated off of an <code>unsigned int</code> describing the number of quadrature samples to be used per triangle when computing integrals.
<UL>
<LI> The mass matrix for extrinsic vector-fields can be obtained by invoking the member function:
<BLOCKQUOTE><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::template massMatrix&lt;QuadratureSamples&gt;( void )</CODE></BLOCKQUOTE>
<LI> The stiffness matrix for extrinsic vector-fields, with respect to the connection energy, can be obtained by invoking the member function:
<BLOCKQUOTE><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::template stiffnessMatrix&lt;QuadratureSamples&gt;( void )</CODE></BLOCKQUOTE>
<LI> The stiffness matrix for extrinsic vector-fields, with respect to a particular component (or combination of components) of the covariant derivative, can be obtained by invoking the member function:
<BLOCKQUOTE><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::template stiffnessMatrix&lt;QuadratureSamples,Component&gt;( void )</CODE></BLOCKQUOTE>
with <code>Component</code> an <code>unsigned int</code> flag indicating whether the divergence, curl, and/or anti-holomorphic components should be used in defining the stiffness. Possible values for <code>Component</code> include:
<UL>
<LI><CODE>SimplexProcessing::StiffnessComponent::Divergence()</CODE>
<LI><CODE>SimplexProcessing::StiffnessComponent::Curl()</CODE>
<LI><CODE>SimplexProcessing::StiffnessComponent::AntiHolomorphic()</CODE>
<LI><CODE>SimplexProcessing::StiffnessComponent::Connection()</CODE>
<LI><CODE>SimplexProcessing::StiffnessComponent::Holomorphic()</CODE>
<LI><CODE>SimplexProcessing::StiffnessComponent::Hodge()</CODE>
<LI><CODE>SimplexProcessing::StiffnessComponent::Killing()</CODE>
</UL>
</UL>

<DT><B>Vector-Field Representation</B>
<DD> While a vector-field can be represented as an <code>std::vector</code> of <code>MishaK::Point&lt;double,3&gt;</code> objects describing the Euclidean coordinates of the (tangent) vector-field values at the vertices, the code also supports a more abstract functional representation. In particular, if an object <code>VectorField</code> satisfies the concept:
<BLOCKQUOTE>
<CODE>
concept HasMeshVectorField = requires( const VectorField f , size_t idx , MishaK::Point&lt;double,2&gt; p ) { { f[idx](p) } -> std::same_as&lt; MishaK::Point&lt;double,3&gt; &gt; };
</CODE>
</BLOCKQUOTE>
behaving as an array of functions on triangles (with triangles indexed by <code>idx</code> and triangle positions described by <code>p</code>) that return the Euclidan representation of a vector, then the dual representation of the vector-field represented by such an object as a vector of size <code>3&middot;|V|</code>, obtained by integrating against the extrinsic Phong-Rodrigues basis vector-fields, can be obtained by invoking the member function:
<BLOCKQUOTE><CODE>
Eigen::VectorXd EmbeddedPhongMesh&lt;2&gt;::template massVector&lt;QuadratureSamples&gt;( VectorField && VF )
</CODE></BLOCKQUOTE>

<DT><B>Example: Computing the Killing Vector-Field</B>
<DD>Given vertex positions, vertex normals, and triangle incidence, the (most) Killing vector-field can be obtained as follows:
<UL>
<LI>Construct an <code>EmbeddedPhongMesh&lt;2&gt;</code> object,
<LI>Compute the extrinsic mass and Killing stiffness energy matrices,
<LI>Restrict those to tangent energy matrices,
<LI>Solve the generalized eigen-value problem to get the tangent vector-field with smallest Killing energy, and
<LI>Prolong to obtain the Euclidean representation of the vector-field at each vertex.
</UL>
Assuming three quadrature samples per triangle, this can be implemented as follows:
<PRE><CODE>static const unsigned int QuadratureSamples = 3;

// Set the triangle mesh information
std::vector&lt; MishaK::Point&lt;double,3&gt; &gt; vertices , normals;
std::vector&lt; MishaK::SimplexIndex&lt;2&gt; &gt; triangles;
// ...

// Initialize the mesh
MishaK::SimplicialMesh::EmbeddedPhongMesh&lt;2&gt; mesh( vertices , normals , triangles );

// Compute the extrinsic/Euclidean system matrices
Eigen::SparseMatrix&lt;double&gt; mass_euclidean = mesh.template massMatrix&lt; QuadratureSamples &gt;();
Eigen::SparseMatrix&lt;double&gt; stiffness_euclidean = mesh.template stiffnessMatrix&lt; QuadratureSamples , MishaK::SimplexProcessing::StiffnessComponent::Killing() &gt;();

// Compute the prolongation matrix from tangent to Euclidean coordinates
Eigen::SparseMatrix&lt;double&gt; prolongation = mesh.tangentProlongation();

// Comptue the tangent system matrices
Eigen::SparseMatrix&lt;double&gt; mass_tangent = prolongation.transpose() * mass_euclidean * prolongation;
Eigen::SparseMatrix&lt;double&gt; stiffness_tangent = prolongation.transpose() * stiffness_euclidean * prolongation;

// Compute the smallest generalized eigen-vector of the pair of systems {stiffness_tangent,mass_tangent}
Eigen::VectorXd killing_tangent;
// ...

// Compute the corresponding extrinsic/Euclidean vector-field representation
Eigen::VectorXd killing_euclidean = prolongation * killing_tangent;
</CODE></PRE>

</DL>
</OL>
</details>

<hr>
<details>
<summary>
<a name="CHANGES"><b>HISTORY OF CHANGES</b></a><br>
</summary>

<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.00/">Version 1.00</a>:
<ul>
<li> initial source code
</ul>

<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.10/">Version 1.10</a>:
<ul>
<li> added support for computing the Lie bracket, independent of the metric and/or connection.
</ul>

<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.15/">Version 1.15</a>:
<ul>
<li> cleaned up code and renamed functions consistently
</ul>

</details>


<hr>
<a href="https://www.cs.jhu.edu/~misha">HOME</a>
