# MATLAB-Voro: MATLAB mex files for the Voro++ software library
MATLAB mex files to run [Voro++](https://math.lbl.gov/voro++/) for calculating [Laguerre diagrams/power diagrams](https://en.wikipedia.org/wiki/Power_diagram) in 2D and 3D periodic and non-periodic domains.

## Installation ##

* Download and change directory to the top-level.
* In MATLAB, run the scripts ``makemex.m`` and ``makemex2d.m``

## Examples ##

See the MATLAB live scripts ``Examples2d.mlx`` and ``Examples3d.mlx`` or the corresponding PDF files ``Examples2d.pdf`` and ``Examples3d.pdf``.

## Voro++ ##

Insert: Cite the Voro++ website, repository and paper. <https://math.lbl.gov/voro++/>, <https://github.com/chr1shr/voro>

## License ##

Insert: Add info on license.

## Software limitations ##

* The domain of the Laguerre tessellation must be a box (a rectangle in 2D or a cuboid in 3D). In particular, the code does not compute unbounded Laguerre tessellations or, e.g., Laguerre tessellations of polygonal domains. 
* 3D code: If the Laguerre tessellation is non-periodic, then all the seeds must lie inside the bounding box, otherwise they are simply discarded. This limitation does not apply to the 3D periodic code or the 2D code.
* 3D code: The Laguerre tessellation must be either non-periodic or triply periodic (periodic in all 3 directions). This limitation does not apply to the 2D code, where you can choose the diagram to be periodic in the x- or y-directions (Laguerre tessellation of a flat cylinder), both directions (Laguerre tessellation of a flat torus) or neither direction.

We plan to address some of these limitations in a future update, as well as adding more sophisticated visualisation tools.

## Related software ##

* Power Diagrams, MATLAB File Exchange: <https://www.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams>
* Fast Bounded Power Diagram, MATLAB File Exchange: <https://www.mathworks.com/matlabcentral/fileexchange/56633-fast-bounded-power-diagram> 
* Geogram: <https://github.com/BrunoLevy/geogram>

## Main contributors ##
This repository was created by

* [Steve Roper](https://www.gla.ac.uk/schools/mathematicsstatistics/staff/stevenroper/#), University of Glasgow
* [David Bourne](http://www.macs.hw.ac.uk/~db92/), Heriot-Watt University and the Maxwell Institute for Mathematical Sciences

This code was used to create Laguerre tessellations for the following papers:

* Egan, C.P., Bourne, D.P., Cotter, C.J., Cullen, M.J.P., Pelloni, B., Roper, S.M. & Wilkinson, M. (2022) A new implementation of the geometric method for solving the Eady slice equations. [Preprint](https://cvgmt.sns.it/paper/5457/)
* Bourne, D.P., Kok, P.J.J., Roper, S.M. & Spanjer, W.D.T. (2020) Laguerre tessellations and polycrystalline microstructures: A fast algorithm for generating grains of given volumes, *Philosophical Magazine*, 100, 2677-2707. <https://doi.org/10.1080/14786435.2020.1790053>
