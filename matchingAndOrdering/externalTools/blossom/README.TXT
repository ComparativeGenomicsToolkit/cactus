BLOSSOM V - implementation of Edmonds' algorithm for computing a minimum cost perfect matching in a graph
Version 2.02
http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html

Details of the implementation are described in

	Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
	In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

Please send comments to vnk@adastral.ucl.ac.uk.
If you use this software for research purposes, you should cite the aforementioned paper in any resulting publication.

##################################################################

License & disclaimer:

    Copyright 2008-2009 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk)

    This software can be used for research and evaluation purposes only. Commercial use is prohibited.
    Public redistribution of the code or its derivatives is prohibited.
    If you use this software for research purposes, you should cite the following paper in any resulting publication:
        Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
        In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

    For commercial use of the software not covered by this agreement, please contact the author.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


##################################################################

COMPILATION:

In unix, type "make". The code should also compile in Windows with Microsoft Visual C++ compiler.
Tested on 32-bit machines.

##################################################################

USAGE:

See PerfectMatching.h for interface functions. Alternatively,
compile the code and run ./blossom5 [options] (see USAGE.TXT for details).

The code also allows solving complete geometric instances; see GEOM/GeomPerfectMatching.h
for interface functions.

##################################################################

PARAMETERS:

For many types of problems, default parameters should be ok.
But for certain problems (such as structured geometric instances)
you may consider setting 
dual_LP_threshold=0.005 (for example). This corresponds to calling
./blossom5 -m0.005
Type PerfectMatching::REAL in PerfectMatcing.h should then be set to double.

##################################################################

EXTERNAL PACKAGE:

When solving complete geometric instances you need to provide the initial subset of edges.
It may be desirable to use Delaunay triangulation for this purpose.
Then you need to download the "Triangle" package of J. R. Shewchuk from 
	http://www.cs.cmu.edu/~quake/triangle.html
and extract it to the directory "triange".
Alternatively, you can use nearest neighbours initialization (this does not require external packages).

