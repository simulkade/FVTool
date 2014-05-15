# FVTool: Finite volume toy toolbox for Matlab
This is a finite volume (toy) toolbox for chemical/petroleum engineers. 
Right now, it can solve a transient convection-diffusion equation with variable velocity field/diffusion coefficients. The discretization schemes 
include:
  * central difference diffusion term
  * central difference convection term
  * upwind convection term
  * TVD convection term
  * transient term
  * Dirichlet, Neumann, Robin, and periodic boundary conditions
## How to start
Download the package, start matlab, and run
   FVToolStartUp
   
You can see a very early document and tutorial by typing (in matlab command line)
   showdemo FVTdemo
   
## Inspiration
I started writing this tool after playing with [FiPy] (http://www.ctcms.nist.gov/fipy/), an amazing python-based finite volume solver. 
This matlab solver is not a clone, and indeed very limited compared to FiPy.
I wrote it to have a very handy tool for testing new ideas (new mathematical models) by solving them in 1D uniform Cartesian grids. 
Then I extended the code to 
  * 1D axisymmetric (radial)
  * 2D Cartesian
  * 3D Cartesian
  * 2D axisymmetric
  
I have overloaded some of the matlab operators to simplify the switch from 1D codes to 2D and 3D.

## Examples
There are a few simple examples in the [Tutorial] (https://github.com/simulkade/FVTool/tree/master/Examples/Tutorial) folder. 
You can also find a few more advanced examples (water injection into a heterogeneous oil field, two nonlinear PDE's, coupled 
fully implicit solution) in the [Advanced] () folder.

## Documents
comming soon

## But Matlab is not a free software?
You can use the code in [octave] (http://www.gnu.org/software/octave/). 
It's not fully compatible because you cannot overload operators for octave builtin types. 
If you avoid using the overloaded operators, it works fine. I've written/am writing the code in [Julia] (http://julialang.org/). It goes well, but 
the visualization and sparse solvers of Julia are not mature yet.

## Other codes
I have a PVT tool (phase equilibrium calculations) in matlab, which is not ready to be shared yet. The code must be commented and revised.
That's the reason that you get this warning durin the start-up.
