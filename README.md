# INCL_python_analysis

Should be relatively simple to use, but this is really really overengineered. 

Essentially just run it like this (if pulled into a folder called INCL_python_analysis):

```bash
python3 INCL_python_analysis -i <input file1> <input file2> etc. -f <additional file1> <additional file2> etc. -t <processing type>
```

this can do these things:

* **kdar:** for a 235 neutrino, gets kdar missing energy comparisons between runs (takes additional files)
* **excitation:** plots just the excitation energy
* **spectral_function:** returns a 2d plot of the spectral function e:p (takes additional files)
* **src:** I cant remember what this does (takes additional files)
* **processing:** returns 8 plots used in validations
* **2p2h:** returns the missing energy for 2p2h events only
* **transparency:** returns comparisons of different transparencies (takes additional files)
* **CCQE_comps:** returns comparisons of excitation energies between samples (takes additional files)

It is definitely too heavy for simple processing, but it is easy to add to once you get the main structure

makes a class which has broken down each event into most of the useful info. this info is then just put into processing scripts, which outputs root files
