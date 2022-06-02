% mexPD_2d.m Help file for mexPD MEX-file
%
%  mexPD_2d.cc - Uses the c++ library voro++ to generate 2D Laguerre diagrams
%  in rectangular domains, with and without periodicity.
%
%  The standard calling syntax is
%
%     [v,t,xc]=mexPD_2d(box,x,w,periodic_x,periodic_y)
%
%   where box is a 1x4 array [xmin ymin xmax ymax] specifying the container size
%         x   is an Nx2 array specifying the locations of the generators
%         w   is an Nx1 array specifying the weights
%         periodic_x is a boolean variable, if 'true' the computed diagram is periodic in the x-direction
%         periodic_y is a boolean variable, if 'true' the computed diagram is periodic in the y-direction
%
%   The output is
%         v     an Nx1 array of cell volumes
%         t     an Nx1 array of cell transport costs
%         xc    an Nx2 array of cell centroids
%
%   If periodic is ommitted the assumption is that the periodic is 'false'
%   If w is ommitted then the assumption is that the weights are all zero
%   If box is ommitted then the assumption is that box=[0 0 1 1] so the
%   domain is the unit square
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
