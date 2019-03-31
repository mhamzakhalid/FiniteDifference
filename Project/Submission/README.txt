

Instructions:

There are runnable 3 files in this project:


testorder2_4:
	Produces solution plot. Can be changed to run different order FD schemes, and plotting their convergence plots, by changing which line is uncommented under %%TASK.


femtest:
	Produces FEM scheme plot.


testMultigrid:
	Uses PCG (preconditioned conjugate gradient) to solve the problem. Under %%TASK, change which line in uncommented to see the convergence for the other methods: CG(regular conjugate gradient) and MGV (multigrid).


In MATLAB, simply click "Run" to have the script return the relevant plots.
The other files included contain dependencies for the 3 runnable files.
