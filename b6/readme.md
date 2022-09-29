## OLD README

This ketchup version b6, submitted as online supplementary material with 
the article 'Vlasov simulations of parallel potential drops' in Annales 
Geophysicae.

The main program is in the file fvlasov.fpp. Subroutines that are specific 
to the model of the L=7 shell used in that article are found in ModelL7.fpp.
With the setup in the makefile included here, the main program executable 
will be l7/ketchup.

It is possible to generate new distribution functions in dump files one can 
start from by using the program in Regeneration.fpp. This is useful for 
refining the grid. The executable is l7/regenerate_ketchup. 

If the distributions are regenerated, it may be necessary to update the 
value of the stored maximum f0 (see article). This can be done with the 
maxf0update program.

The directories l7/testcase1 and l7/testcase2 contain the input files that 
are needed to run the program. The first case can be run with ketchup, 
and the second can then be started from the dump files produced by 
the first by running regenerate_ketchup before running ketchup itself. 
If you do not want to load initial densities, temperatures, etc from files, 
the files named *.inp can be removed, as in l7/testcase3, which starts from 
an empty state. 

Matlab m-files for file conversion and plotting are found in directory 
m_files.
