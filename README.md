
# Graph Partitioning
=====================================

This code implements the algorithms presented in the paper "Solving Graph Partitioning on Sparse Graphs: Cuts, Projections, and Extended Formulations" by Demetrios V. Papazaharias and Jose L. Walteros. The code, which is implemented in C++, includes a variety of formulations and cutting-plane subroutines applied to solve graph partitioning problems. The implementation uses the libraries of commercial solver Gurobi. 


Usage
---------

### Compiling the code
The library provides a simple ``Makefile`` in the root folder to compile and run the experiments described in the aforementioned paper. The library has been tested on Linux (Ubuntu, Manjaro, and CentOS) and macOS. For the general compilation executing the basic Linux command ``$make`` in the root folder will suffice. Importantly, the user requires to update the variable ``GRBPATH`` to reflect the path of the folder containing the Gurobi library, as well as the corresponding flags of the version installed in the computer. 

####[For the MPC reviewers. The library is already compiled in the VM, so they can omit the compilation step. In any case, the ``GRBPATH`` and flags are also properly set, so if the reviewers want to recompile the code they can also do so. In such a case, we suggest running ``$make clean`` first].

### Input file format
The code can receive as input graphs given in the edge list format. See the dat folder for several examples. In fact all the data files with the instances used in the paper above are included in that folder.

### Running the code
To replicate the experiments of the paper we have included multiple Bash Shell Scripts that call the proposed algorithms and provide the required input parameters for each run. Importantly, the paper compiles the computational results in 9 tables (Table 2, ..., Table 10). To reproduce the experiments we used to create each table, the user can run any of the 7 scripts available in the ``./run`` folder (i.e., ``table2.sh, table3.sh, table4.sh, table5-6.sh, table7-8.sh, table9.sh, and table10.sh``). The best and suggested option to run such scripts is via the ``Makefile``. This option is preferred as it guarantees that the path of the input files is set correctly. To run the experiments associated with each table, the user can run the command:

	$make table<number_of_the_table>
For example to run the experiments for table 2 use:

	$make table2

Importantly, tables 5 and 6 use the same experiments, so as tables 7 and 8. Thus, the commands to run the experiments for those tables are

	$make table5-6
	$make table7-8

Finally, if the user wants to run all the experiments directly, they can use the command

	$make run_all

This will call all the aforementioned commands directly. 

Individual runs that call the algorithms can also be ran. The following is the input information. Depending on formulation selected, addition arguments may be required.

	argv[1]: filetype - input file type edge list representation (-e) or adjacency list representation (-a)
	argv[2]: weighted - unweighted (-u), edge weighted (-ew), randomly generated weights (-w)
	argv[3]: filename - the instance we wish to solve"
	argv[4]: formulation - Options : -TRI, -TRILP, -FLOW, -FLOW+, -PATH, -PATH+, -TCFLP, -TDPLP, -TDPBLP


### Output files
If any of the given scripts is used, the corresponding output file will appear in the folder ``./out``. The output of each run may be different, thus the first line(s) of each of each output files contains a header with the output information. Averages for the tables can be calculated using the raw output files. 

Contact information
--------------------
For inquiries about the implementation and results, please contact:

####Jose L. Walteros, Ph.D.  
Assistant Professor  
Department of Industrial and Systems Engineering  
University at Buffalo, The State University of New York  
413 Bell Hall, Buffalo, NY 14260  
Tel: 716-645-8876, Fax:716-645-3302  
Email: josewalt@buffalo.edu 


Terms and conditions
--------------------

Please feel free to use these codes. We only ask that you cite: 

####"Solving Graph Partitioning on Sparse Graphs: Cuts, Projections, and Extended Formulations" by Demetrios V. Papazaharias and Jose L. Walteros

_This library was implemented for academic purposes. You may use the library with no limitation for academic, commercial, and/or leisure applications without the need an explicit permission, other than the one granted by the MIT type license included in the library. We kindly ask you to aknowledge the library by citing the above paper in any publication that uses this resources._

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.






