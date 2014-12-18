AFRODITE
========

A GEANT4 simulation of the AFRODITE vault at iThemba LABS, composed of the AFRODITE array of HPGe clover detectors along with various ancillary detectors

A CAD model import interface called CADMesh (authored primarily by Christopher Poole) is used within this simulation for implementing complex geometries. This needs to be obtained from the following GitHub repository:
https://github.com/christopherpoole/CADMesh

Alternatively, one could comment out the relevant CADMesh associated code and use only the hard-coded geometrical objects. It should be noted that an effort has been made to hard-code the geometries in GEANT4 when possible as this has computational advantages.