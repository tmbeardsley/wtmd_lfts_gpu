# Well-tempered metadynamics applied to Langevin field-theoretic simulations of diblock copolymers on GPUs

See https://www.tbeardsley.com/projects/lfts/wtmd for a detailed discussion of this project.<br>

## Required Dependencies:
GSL - GNU Scientific Library (https://www.gnu.org/software/gsl/)<br>
CUDA Toolkit (https://developer.nvidia.com/cuda-toolkit/)<br>

## Compiling:
Two methods of compiling the program are available:<br>
<ol>
  <li><b>comp.sh</b>
    <br>
    A simple bash script to create a 'build' directory containing the compiled program code: wtmd-lfts-gpu.<br><br>
    On a Linux system, run the bash script from the top directory via:<br>
    <b>sh comp.sh</b>
    <br><br>
  </li>
  <li><b>CMake</b>
    <br>
    CMakeLists.txt specifies the required commands for CMake to create (and run) Makefiles, which create a 'build' directory and compile the program code as: wtmd-lfts-gpu.<br><br>
    From the top directory, run: <br>
    <b>cmake -B build</b><br>
    <b>cmake --build build</b>
  </li>
</ol>

## Running the program
After compilation the executable file, wtmd-lfts-gpu, resides in the 'build' directory. Two input files are required and must be supplied to the executable at the command line. The first specifies the standard parameters and input fields for a Langevin field-theoretic simulation (L-FTS). The second specifies parameters and input potentials related specifically to performing well-tempered metadynamics (WTMD). Examples of these two files are contained in the 'input_files' folder. 
For example, from the top level of the directory tree the program could be run via: <br><br>
<b>./build/wtmd-lfts-gpu ./input_files/input ./input_files/bias_in</b>

## Input Files
The input_files directory contains example input files that can be supplied to the program from the command line.

### L-FTS Input file

#### File Format
Line 1: <em>N NA XN C Ndt isXeN</em><br>
Line 2: <em>mx my mz Lx Ly Lz</em><br>
Line 3: <em>n_eq n_st n_smpl save_freq loadType</em><br>
Lines 4->(M+3): W-(r)<br>
Lines (M+4)->(2M+3): w+(r)<br>

Note: A real-space position r = (x,y,z) corresponds to a mesh point position r_m = (i,j,k), where i=0->mx-1, j=0->my-1 and k=0->mz-1 are integers. The elements of the fields, W-(r) and w+(r), are then written in ascending order of the row-major index: p = mx\*(i\*my+j)+k.

#### Parameter Descriptions
<em>N</em> is the number of monomers in a single polymer chain (integer).<br>
<em>NA</em> is the number of monomers in the A-block of a polymer chain (integer).<br>
<em>XN</em> is the interaction strength between A anD B-type monomers (double).<br>
<em>C</em> is the square root of the invariant polymerisation index, Nbar (double).<br>
<em>Ndt</em> is the size of the time step in the Langevin update of W-(r) (double).<br>
<em>isXeN</em> instructs the program whether the parameter XN is in terms of bare (isXeN=0) or effective (isXeN=1) chi (integer).<br>
<em>mx, my, mz</em> are the number of mesh points in the x, y, and z dimensions of the simulation box (integers).<br>
<em>Lx, Ly, Lz</em> are the dimensions of the simulation box (in units of the polymer end-to-end length, R0) in the x, y, and z dimensions (doubles).<br>
<em>n_eq</em> is the number of langevin steps performed to equilibrate the system (integer).<br>
<em>n_st</em> is the number of langevin steps performed after equilibration has ended, during which statistics are sampled (integer).<br>
<em>n_smpl</em> is the number of steps between samples being taken in the statistics period (integer).<br>
<em>save_freq</em> is the number of steps between saving outputs to file.<br>
<em>loadType</em> instructs the program whether to load the W-(r) and w+(r) fields from the proceeding file lines (loadType=1), start from a disordered state (loadType=0) or start from a (300) lamellar phase (loadType=2).<br><br>
M = (mx\*my\*mz) is the total number of mesh points, such that the proceeding 2*M lines of the file can hold W-(r) and w+(r) fields to load.<br>

### Bias Field Input File

#### File Format
Line 1: <em>ell kc sigma_Psi DT Psi_min mPsi dPsi update_freq read_bias</em><br>
Lines 2->(mPsi+1): Psi u(Psi) up(Psi) I0(Psi) I1(Psi)<br>

#### Parameter Descriptions
<em>ell</em> is the constant, l, used in the definition of the order parameter (double).<br>
<em>kc</em> is the wavevector cutoff (constant) in order parameter (double).<br>
<em>sigma_Psi</em> is the width of Gaussians added to the bias potential (double).<br>
<em>DT</em> controls the rate at which the amplitude of Guassians is reduced (double).<br>
<em>Psi_min</em> is the lowest value of Psi for which the bias potential is collected (double).<br>
<em>mPsi</em> is the total number of mesh points covering the range of Psi being investigated (integer).<br>
<em>dPsi</em> is the distance between adjacent mesh points in Psi (double).<br>
<em>update_freq</em> is the number of Langevin steps between updates of the bias potential (integer).<br>
<em>read_bias</em> is a flag indicating whether an initial bias potential should be read from file (read_bias=1) or not (read_bias=0) (integer).
<br><br>
Lines 2->(mPsi+1) are only read if read_bias=1, otherwise the bias potential starts with all elements equal to zero.<br>

## Output Files
#### w_eq_<step_number>
The state of the W-(r) and w+(r) fields at simulation step number <step_number> during the equilibration period. First three lines are simulation parameters so it can be used as an input file.<br>

#### w_st_<step_number>
The state of the W-(r) and w+(r) fields at simulation step number <step_number> during the statistics gathering period. First three lines are simulation parameters so it can be used as an input file.<br>

#### phi_eq_<step_number>
The state of the phi-(r) and phi+(r) fields at simulation step number <step_number> during the equilibration period.<br>

#### phi_eq_<step_number>
The state of the phi-(r) and phi+(r) fields at simulation step number <step_number> during the statistics gathering period.<br>

#### bias_st_<step_number>
Output file in the same format as the WTMD input file, containing the current state of the bias potential after step_number Langevin steps in the statistics gathering phase.
