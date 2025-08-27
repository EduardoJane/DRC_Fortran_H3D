# DRC_Fortran_H3D
This is an extension to the open-source [Delft Research Controller (DRC)](https://github.com/TUDelft-DataDrivenControl/DRC_Fortran) (see [Referencing](#referencing) for more information) prepared to have full compatibility with the external controller module for Actuator Lines in the High-Order CFD solver [HORSES3D](https://github.com/loganoz/horses3d). The most relevant additions include the ability to write the controller state to disk (filter and PID cumulative terms among other state variables) and to accordingly read the written state file at restart. Both operations are triggered by the custom controller status flags defined in the HORSES3D external controller module. 

To compile the external library (*.so) controller file, simply use `make` with the provided [Makefile](./Source/makefile), making sure the compiler and its flags are the same used to compile HORSES3D to ensure compatibility.

The original README for the DRC repository can be found below.

# DRC_Fortran wind turbine baseline controller (original README)
Delft Research Controller (DRC) baseline wind turbine controller, using the Bladed-style DISCON interface used by, e.g., OpenFAST, Bladed (versions 4.5 or earlier) and HAWC2.

## Introduction
The Delft Research Controller (DRC) provides an open, modular and fully adaptable baseline wind turbine controller to the scientific community. New control implementations can be added to the existing baseline controller, and in this way, convenient assessments of the proposed algorithms is possible. Because of the open character and modular set-up, scientists are able to collaborate and contribute in making continuous improvements to the code. The DRC is being developed in Fortran and uses the Bladed-style DISCON controller interface. The compiled controller is configured by a single control settings parameter file, and can work with any wind turbine model and simulation software using the DISCON interface. Baseline parameter files are supplied for the NREL 5-MW and DTU 10-MW reference wind turbines.

## Using the DRC for Bladed
If you want to use the controller with DNV GL Bladed v4.5 or earlier (which still has support for the DISCON external controller interface), do the following:
1. Be sure to use and place the 32-bit DLL in the same folder as where you put your project .$PJ-file
2. Copy in that same folder the DISCON.IN controller configuration file
3. Set-up the 32-bit DLL as an external controller (Control -> Discrete External Controller -> Define...)
3. Open the DISCON.IN file with a text editor and copy its entire contents in the "External controller data:" section (Control -> Discrete External Controller -> Define...)
4. Run a "Power Production" simulation

## Referencing
When you use the DRC in any publication, please cite the following paper:
* Mulders, S.P. and van Wingerden, J.W. "Delft Research Controller: an open-source and community-driven wind turbine baseline controller." Journal of Physics: Conference Series. Vol. 1037. No. 3. IOP Publishing, 2018. [Link to the paper](https://iopscience.iop.org/article/10.1088/1742-6596/1037/3/032009/meta)