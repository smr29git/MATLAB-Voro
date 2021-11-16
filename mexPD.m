% mexPD.m Help file for mexPD MEX-file
%  mexPD.cc - Uses the C++ library voro++ to generate 3D Laguerre diagrams
%  in cuboidal domains, with and without periodicity.
%
%  The standard calling syntax is
%
%     [v,t,xc]=mexPD(box,x,w,periodic)
%
%   where box is a 1x3 array [l1 l2 l3] specifying the container size
%         x   is an Nx3 array specifying the locations of the generators
%         w   is an Nx1 array specifying the weights
%         periodic is a boolean variable, if 'true' the computed diagram is
%         periodic if 'false' the computed diagram is not periodic
%
%   The output is
%         v     an Nx1 array of cell volumes
%         t     an Nx1 array of cell transport costs 
%         xc    an Nx3 array of cell centroids
%
%   If periodic is ommitted the assumption is that the periodic is 'false'
%   If w is ommitted then the assumption is that the weights are all zero
%   If box is ommitted then the assumption is that box=[1 1 1] so the
%   domain is the unit cube
%
%   It is therefore possible to call mexPD with the following arguments
%
%     [v,t,xc]=mexPD(x);
%
%     [v,t,xc]=mexPD(box,x);
%     [v,t,xc]=mexPD(x,w);
%     [v,t,xc]=mexPD(x,periodic);
%
%     [v,t,xc]=mexPD(box,x,w);
%     [v,t,xc]=mexPD(box,x,periodic);
%     [v,t,xc]=mexPD(x,w,periodic);
