# FATTYRIOT Top-level function for 2012 ISMRM Fat-Water Challenge

Requirements:
=============
* Mac OS X (tested with Lion) or Windows (tested with 32-bit WinXP) or i386 LINUX (tested on Ubuntu 8.04 LTS)
* Tested with MATLAB R0210a, MATLAB 2012a and MATLAB 2013a

Installation:
=============
* Run `setup_FattyRiot.m` keeping it in the same folder as `FattyRiot.m`
  - Adds the parent folder of `FattyRiot.m` to the MATLAB path
  - All other needed paths are added dynamically by `FattyRiot.m`
* LINUX users need to unzip the file `./FattyRiot_toolbox/berglund/QPBO/DixonApp/LINUX/DixonApp_LINUX.exe.zip` (GITHUB does not allow files more than 100 MB in size)

Usage:
======
* `FattyRiot.m` calls the FattyRiot algorithm using syntax

`FW = FattyRiot(imDataParams);`

  Input: ISMRM fat-water toolbox structure
  - imDataParams.images                : acquired images, array of size [nx, ny, nz, ncoils, nTE]
  - imDataParams.TE                    : echo times (in seconds), vector of length nTE
  - imDataParams.FieldStrength         : (in Tesla), (default 1.5)
  - imDataParams.PrecessionIsClockwise : ==1 is clockwise, ~=1 is counterclockwise, (default 1) 
  - imDataParams.mask                  : logical [nx, ny, nz], ( default true([nx, ny, nz]) )

  Output:
  - FW (fat and water magnitude images) size [nx, ny, 2*nz]
  
Example:
========
* `test_FattyRiot.m`
  - Tests FattyRiot on all 10 Phase I and 7 Phase II cases from the 2012 ISMRM Challenge
  - should be run from the same folder as `FattyRiot.m`
