# E.coli_wholecellmodel
This is an analytical plotting tool to assist whole-cell models of E. coli. It needs to be used on a supercomputer with  E. coli whole-cell model. Download Runningfile to the E. coli whole cell model  folder you installed.

Running:
python run.py path1 path2 path3…… exportpath
each path is separated by a space. 
For example: python run.py \ /newhome/wo20754/wholecell3/wcEcoli/out/20210523.151321__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/wildtype_000000/000000/generation_000000/000000/simOut/ \ /newhome/wo20754/wholecell3/wcEcoli/out/20210512.154744__geneknockdown_Run_2gen3seedVARIANT=geneKnockdown/wildtype_000007/000000/generation_000000/000000/simOut/  \
/newhome/wo20754/wholecell3/wcEcoli/out/project/exportpath/
Then you will get the figure from the output of path1 and path2 which has the same X range and Y range.
