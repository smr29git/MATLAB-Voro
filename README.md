# MATLAB-Voro: MATLAB mex files for the Voro++ software library
MATLAB mex files to run [Voro++](https://math.lbl.gov/voro++/) for calculating [Laguerre diagrams/power diagrams](https://en.wikipedia.org/wiki/Power_diagram) in 2D and 3D periodic and non-periodic domains.

## Installation ##

* Download and change directory to the top-level.
* In MATLAB, run the script ``makemex.m``

## Examples ##

Live scripts & PDF

## Software limitations ##

* The domain of the Laguerre tessellation must be a box (a rectangle in 2D or a cuboid in 3D). In particular, it does not compute unbounded Laguerre tessellations or, e.g., Laguerre tessellations of polygonal domains. 
* 3D code: If the Laguerre tessellation is non-periodic, then all the seeds must lie inside the bounding box, otherwise they are simply discarded. This limitation does not apply to the 3D periodic code or the 2D code.
* 3D code: The Laguerre tessellation must be either non-periodic or triply periodic (periodic in all 3 directions). This limitation does not apply to the 2D code, where you can choose the diagram to be periodic in the x- or y-direction (Laguerre tessellation of a flat cylinder), both directions (Laguerre tessellation of a flat torus) or neither direction.

We plan to address some of these limitations in a future update, as well as adding more sophisticated visualisation tools.


## Related software ##

* Voro++: <https://math.lbl.gov/voro++/>, <https://github.com/chr1shr/voro>
* Power Diagrams, MATLAB File Exchange: <https://www.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams>
* Fast Bounded Power Diagram, MATLAB File Exchange: <https://www.mathworks.com/matlabcentral/fileexchange/56633-fast-bounded-power-diagram> 
* Geogram: <https://github.com/BrunoLevy/geogram>

## Things to add before publishing ##
* License file and info
* Cite the Voro++ website, repository and paper
* Add the authors of this page (Steve and David) and cite any relevant papers where we've used this code (just the Tata paper for now, Mason's future papers)
