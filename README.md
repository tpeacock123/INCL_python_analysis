# INCL_python_analysis

Should be relatively simple to use, but this is really really overengineered. 

Essentially just run it like this (if pulled into a folder called INCL_python_analysis):

python3 INCL_python_analysis -i < input file1 > < input file2 > etc. -f < additional file1 >  < additional file2> etc. -t 
< processing type >


this can do these things:
kdar:              for a 235 neutrino, gets kdar missing energy comparisons between runs (takes additional files)\n"
excitation:        plots just the excitation energy\n"
spectral_function: returns a 2d plot of the spectral function e:p (takes additional files)\n"
src:               I cant remember what this does (takes additional files)\n"
processing:        returns 8 plots used in validations\n"
2p2h:              returns the missing energy for 2p2h events only\n"
transparency:      returns comparisons of different transparencies (takes additional files)\n"
CCQE_comps:        returns comparisons of excitation energies between samples (takes additional files)

It is definitely too heavy for simple processing, but it is easy to add to once you get the main structure

makes a class which has broken down each event into most of the useful info. this info is then just put into processing scripts, which outputs root files