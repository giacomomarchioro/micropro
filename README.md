# micropro
Python module for the analysis of microprofilometry data. This module was produced during my PhD thesis and allows to open and analyze the data produced using the microprofilometric setup implemented at OpDaTeCH laboratory. 

# Installation 
  ```
  pip install git+https://github.com/giacomomarchioro/micropro
  ```
# Basic usage
  
  ```python
  from micropro import ns
  mysurface = ns.load(path) # Create an instance for the surface
  mysurface.mfilter() # filter bad values
  mysurface.subtractplane(order=1) # fit plane and subtract it form the surface
  mysurface.plot() # plot the results
  ```


This work was supported by the Scan4Reco project funded by EU Horizon 2020 Framework Programme for Research and Innovation under Grant Agreement no 665091.

