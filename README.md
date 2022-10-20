# ketchup

This is a time-dependent electrostatic Vlasov simulations of double layers (DL) in an auroral flux tube. 

The code is seperated in two folders:
- `b6/` which is almost identical to the code published as online supplementary material with the article of [Gunell et al. (2013)](https://doi.org/10.5194/angeo-31-1227-2013).
- `MI_coupling/` which is a modified version of the `b6/` code that supports an input from a file at the ionospheric boundary. The ionospheric response can typically be calculated using the time-dependent electron transport [AURORA](https://github.com/egavazzi/AURORA) model described in [Gustavsson, B. (2022)](https://doi.org/10.1029/2019JA027608)

<br />

## Installation
Download and extract the .zip or clone the repository (e.g. using Git).

You should also make sure that the [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.derk94) is installed on your machine.

<br />

## Get started
To be able to compile and run the Fortran code from a terminal, you need to export some environment variables. The necessary commands are gathered in the shell script `ketchup/init_bash.sh` so that you just need to run it in your terminal.


The best is to start with the `b6/` code. To compile the Fortran code (all the *.fpp* files) situated there, use the Makefile with:
```
$> make all
```
This will create three executable: 
- `b6/l7/ketchup`
- `b6/l7/regenerate_ketchup`
- `b6/maxf0update`

Then, move into the `l7/` folder, where all experiments are to be run from. There is already a `model` folder situated there, with the following folder structure (required for all experiments):
```
model/
├── dumps/
│   
├── inputb6.m
│ 
├── outp/
│   └── datfiles/...
│
├── poleslevel04s.dat
├── poleslevel04s.m
│
└── regen_par.m
```
- `dumps/ ` is the folder where the program periodically dumps some informations to be used by itself. 
- `inputb6.m` is the file where the experiment parameters are set.
- `outp/...` is the folder where the results are saved.
- `poleslevel04s.dat` and `poleslevel04s.m` are the files defining the coordinate transformations to be used by ketchup.
- `regen_par.m` is a file to be used when starting an experiment from a previous experiment final state.

<br />

For an example, let's' copy the `model/` folder into a `testcase1/` and a `testcase2/` folder.
Move into `testcase1/` and in the `inputb6.m` file, change *startfromdumpfile* to **no**
. Then, you can run the experiment from the terminal using the command:
```
$> mpiexec -ppn X ../ketchup
```
where X is the number of processors threads you want to use (the program will parallelise accordingly). To analyse the data produced by the simulation, you can open Matlab from the `testcase1/` folder and run the `ketchupb6conv.m` script which will convert all the .dat file containing -ASCII data into .mat files. You can then use the `ketchupb6plot.m` file to plot some of the results.

Now let say you want to start another simulation based on the final state of `testcase1/` but with different parameters. To do this, move to `testcase2/`, and set the parameters you want in `inputb6.m`. Then, edit `regen_par.m` and change *FromDir* to **../testcase1/**. Change also *N_procs_new* to the number of processors you will want to use when running the `testcase2` simulation. Once this is done, use the command:
```
$> ../regenerate_ketchup
```
and the simulation is now ready to be run with the same command as before! Just use the same number of processors X than the one you wrote in `regen_par.m`
