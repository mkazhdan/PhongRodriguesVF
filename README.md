<center><h2>Phong-Rodrigues Extrinsic Vector-Field Processing (Version 1.10)</h2></center>
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
<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.10/PRVF.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/Version1.10/PRVF.Source.zip">ZIP</a> <a href="https://github.com/mkazhdan/PhongRodriguesVF">GitHub</a><br>
<B>Data:</B>
<A HREF="https://www.cs.jhu.edu/~misha/Code/PhongRodriguesVF/PRVF.Data.zip">ZIP</A><br>
<b>Older Versions:</b>
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

<DT><B>System Matrix Representation</B></DT>
<DD> After being initialized with the vertices, normals, and triangles of a mesh, different system matrices can be computed.
<UL>
<LI> The method<BR>
<CENTER><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::mass&lt;QuadratureSamples&gt;(void)</CODE></CENTER><BR>
return the symmetric <code>3&middot;|V|&times;3&middot;|V|</code> mass matrix, with <CODE>QuadratureSamples</CODE> the numer of quadrature points per triangle.
<LI> The method:
<CENTER><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::stiffness&lt;QuadratureSamples&gt;(void)</CODE></CENTER><BR>
return the symmetric <code>3&middot;|V|&times;3&middot;|V|</code> stiffness matrix defined by the connection Laplacian.
<LI> Symmetric <code>3&middot;|V|&times;3&middot;|V|</code> matrices defining the stiffness with respect to individual components of the covariant matrix can be obtained invoking the method:<BR>
<CENTER><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::stiffness&lt;QuadratureSamples,CovComponent&gt;(void)</CODE></CENTER><BR>
with <CODE>CovComponent</CODE> an <CODE>enum</CODE> of type <CODE>CovariantComponent</CODE> describing the component of the covariant derivative used to define the stiffness.
<LI> The method:
<CENTER><CODE>Eigen::SparseMatrix&lt;double&gt; EmbeddedPhongMesh&lt;2&gt;::tangentProlongation(void)</CODE></CENTER><BR>
returns the sparse <code>3&middot;|V|&times;2&middot;|V|</code> matrix mapping two degrees of freedom at a vertex, representing the linear combinations of a tangent vector with respect to a chosen frame, to its coordinates in 3D. The frame is chosen to be orthonormal, so that if <CODE>P</CODE> is the matrix returned by the method, the matrix <CODE>P.transpose()*P</CODE> is the <code>3&middot;|V|&times;3&middot;|V|</code> matrix describing the operation of projecting out the normal component.
</UL>
In general, the system matrices should be of size <code>2&middot;|V|&times;2&middot;|V|</code> corresponding to optimizing over tangent vectors. For example, if <CODE>M</CODE> is the mass matrix and <CODE>P</CODE> prolongation, the associated mass matrix defined over tangent vectors will be <CODE>P.transpose() * M * P</CODE>.

<DT><B>System Vector Reprsentation</B></DT>
<DD>
Extrinsic vector fields are represented by <CODE>Eigen::VectorXd</CODE>s of size <code>3&middot;|V|</CODE>, with the <I>x</I>-, <I>y</I>-, and <I>z</I>-coordinates of the <code>v</code>-th vertex found at indices <CODE>3&middot;v</CODE>, <CODE>3&middot;v+1</CODE>, and <CODE>3&middot;v+2</CODE>, respectively.

<DT><B>Abstract Vector-Field Representation</B></DT>
<DD>An abstract vector-field over the mesh is an object of type <CODE>MeshField</CODE> that acts as an array of per-triangle vector-fields of abstract type <CODE>TriangleField</CODE>. Specifically:
<UL>
<LI> The type <CODE>MeshField</CODE> supports a method:<BR>
<CENTER><CODE>TriangleField MeshField::operator[](size_t)</CODE></CENTER><BR>
 giving the restriction of the vector-field to the indexed triangle.
<LI> The type <CODE>TriangleField</CODE> supports a method<BR>
<CENTER> <CODE>Point&lt;double,3&gt; TriangleField(Point&lt;double,2&gt;)</CODE></CENTER><BR>
returning the value of the vector-field within the triangle.
</UL>
One can obtain <CODE>Eigen::VectorXd</CODE> of size <code>3&middot;|V|</CODE> giving the integral of the inner-product of the vector-field with each of the Phong-Rodrigues basis vector-fields by invoking the method:<BR>
<CENTER><CODE>Eigen::VectorXd EmbeddedPhongMesh<2>::dual&lt;QuadratureSamples&gt;(MeshField &amp;&amp;) const</CODE></CENTER>

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

</details>


<hr>
<a href="https://www.cs.jhu.edu/~misha">HOME</a>
