OpenFOAM Project: 

Biofilm Growth Simulation in Streams

Overview:
    
This project simulates biofilm growth in streams using OpenFOAM, a high-fidelity Computational Fluid Dynamics (CFD) software. The study focuses on investigating the impact of various nutrient concentrations and flow rates on biofilm development under different environmental conditions. The results aim to provide insights into biofilm dynamics and help improve water quality management strategies.

Project Structure:

  0/: Initial conditions for the simulation.
  
  constant/: Contains physical properties, mesh files, and boundary conditions.
  
    polyMesh/: The mesh directory containing files defining the computational grid.
    transportProperties: Defines the properties of the fluid.
    biofilmProperties: Custom file to define biofilm-specific properties.
  system/: Configuration files for controlling the simulation.
  
    controlDict: Governs the simulation time, write intervals, etc.
    fvSchemes: Specifies the discretization schemes for the governing equations.
    fvSolution: Contains the solvers and algorithms used for the numerical solution.
    postProcessing/: Scripts and files for post-processing simulation data.
    scripts/: Additional scripts for pre- and post-processing, such as mesh generation and data analysis.
Requirements:

    OpenFOAM v7
    ParaView (for post-processing visualization)
