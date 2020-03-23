*******************************************
Authors
	@ Abhinav Bollam  - 150102075
	@ Saurabh Dhal    - 150102060
	@ Deepak Meghwal  - 150108009 

under the guidance of Dr. Gaurav Trivedi
********************************************

Problem (12c)

Objective: 

	Implementation of Distributed Parallel Cooperative Coevolution based Multi-objective
	optimization algorithm for large scale optimization problems:

Reference paper:

 	https://ieeexplore.ieee.org/abstract/document/7867084

How to run :
	1. Compile main cpp file 
		g++ main.cpp
	2. ./a.out
		
		a) Input for DTLZ type to be given.
			 value from 1-7
		b) Enter Number of variable.   
			Any positive integer value can be given as input. But larger variable values ttake too much time for computation.
		c) Enter number of objective functions.
		d) Enter population size. It can be a positive integer below 1000.

	3. the program runs and outputs run time, number of generations required. (evaluations)
	4. The variable values that are solutions for pareto front solutions are printed to
	 	Var_DTLZ*_*.txt
	5. 	The functional values that are the pareto front values are printed to Fun_DTLZ*_*.txt
	6. The functional values can be plotted using plotfile2d.py for 2 objective functions and 
		plotfile3d.py for 3 objective functions.

All the tabulated results in report are present in variable results and function results. The pictures for 3d are in Result piuctures folder.


