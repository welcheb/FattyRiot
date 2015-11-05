FattyRiot
=========
[![DOI](https://zenodo.org/badge/11902/welcheb/FattyRiot.svg)](http://dx.doi.org/10.5281/zenodo.16741)

FattyRiot is an algorithm for separation of fat and water magnetic resonance images.

This version of FattyRiot is the winning entry (29-Mar-2013 submission with score 9931) of the [International Society for Magnetic Resonance in Medicine (ISMRM) 2012 Challenge on Water-Fat Reconstruction](http://www.ismrm.org/challenge/node/18).

FattyRiot includes algorithms available in the [ISMRM Fat-Water Toolbox v1](http://ismrm.org/workshops/FatWater12/data.htm). 

FattyRiot team members include:

* E Brian Welch (captain), Vanderbilt University
* David S Smith, Vanderbilt University
* Malcolm J Avison, Vanderbilt University
* Johan Berglund, Uppsala University
* Joel Kullberg, Uppsala University
* Håkan Ahlström, Uppsala University

### LICENSE
* FattyRiot is available under a [Creative Commons Attribution-NonCommercial 4.0 International Public License](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode). 
* FattyRiot contains source code and object code from multiple sources.
* FattyRiot is available strictly for non-commerical research purposes.
* See [LICENSE.md](./LICENSE.md) for more details.

### COPYRIGHT
* FattyRiot contains copyrighted material from multiple sources. 
* See [COPYRIGHT.md](./COPYRIGHT.md) for more details.

### SETUP
* See [SETUP.md](./SETUP.md) for more details.

### REFERENCES
* FattyRiot contains algorithms and scientific concepts from many published works. 
* See [REFERENCES.md](./REFERENCES.md) for more details.

### SEE ALSO
* FattyRiot poster presented at 2013 ISMRM: [`2013_ISMRM_(Smith,Welch)_ISMRM_Fat-Water_Challenge_Poster.pdf`](./2013_ISMRM_\(Smith,Welch\)_ISMRM_Fat-Water_Challenge_Poster.pdf)
* GITHUB repository for [fw_i3cm1i_3pluspoint_berglund_QPBO](https://github.com/welcheb/fw_i3cm1i_3pluspoint_berglund_QPBO)

### KNOWN ISSUES
* A bug in the `DixonApp.m` wrapper causes incorrect operation when the water peak location is not entered as 0 ppm. Since this repo is supposed to be a "frozen" version of the highest scoring FattyRiot entry, the bug remains. The bug is corrected in the [fw_i3cm1i_3pluspoint_berglund_QPBO](https://github.com/welcheb/fw_i3cm1i_3pluspoint_berglund_QPBO) repo.  

<hr>
[![Analytics](https://ga-beacon.appspot.com/UA-54485519-2/FattyRiot/README.md)](https://github.com/welcheb/FattyRiot)
