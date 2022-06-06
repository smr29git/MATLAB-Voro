% mexPDall_2d.m Help file for mexPDall MEX-file
%
%
%  mexPDall_2d.cc - Uses the c++ library voro++ to generate 2D Laguerre diagrams
%  in rectangular domains, with and without periodicity.
%
%  The standard calling syntax is
%
%     [area,tr,xc,vfn]=mexPDall_2d(box,x,w,periodic_x,periodic_y)
%
%   where box is a 1x4 array [xmin ymin xmax ymax] specifying the box dimensions
%         x      is an Nx2 array specifying the locations of the generators
%         w      is an Nx1 array specifying the weights
%         periodic_x  is a boolean variable, if 'true' the computed diagram is periodic in the x-direction, if 'false' then it is not
%         periodic_y  is a boolean variable, if 'true' the computed diagram is periodic in the y-direction, if 'false' then it is not
%
%
%   The output is
%         area  an Nx1 array of cell areas
%         tr    an Nx1 array of cell transport costs 
%         xc    an Nx2 array of cell centroids
%         vfn   an Nx2 cell array of cell vertex coordinates, neighbours and remapped seed locations
%               provided the ith cell is non-empty then
%               vfn{i,1} contains an Nv x 2 array of vertex coordinates for the ith cell (Nv number of vertices)
%               vfn{i,2} contains an Nv x 1 array of neighbour indices for the ith cell

%   
%   If periodic_x, periodic_y are omitted the assumption is that the periodic flags are 'false'
%   If w is ommitted then the assumption is that the weights are all zero
%   If box is ommitted then the assumption is that box=[0 0 1 1] so the
%   domain is the unit square
%
%   It is currently possible to call mexPDall_2d with the following arguments
%
%    With one argument
%
%     [area,tr,xc,vfn]=mexPDall_2d(x);
%
%    With two arguments
%
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x);
%     [area,tr,xc,vfn]=mexPDall_2d(x,w);
%
%    With three arguments, per_x=per_y=per when per is specified
%
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,w);
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,periodic);
%     [area,tr,xc,vfn]=mexPDall_2d(x,w,periodic);
%
%    With four arguments, periodic_x=periodic_y
%
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,w,periodic);
%
%    With five arguments
%
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,w,periodic_x,periodic_y);
%