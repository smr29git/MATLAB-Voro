% mexPDall_2d.m Help file for mexPDall MEX-file
%
% Last modified 2021/01/27 SR
%
%  mexPDall_2d.cc - Uses the c++ library voro++ to generate 2D Laguerre diagrams
%  in rectangular domains, with and without periodicity.
%
%  The standard calling syntax is
%
%     [area,tr,xc,vfn]=mexPDall(box,x,w,per_x,per_y)
%
%   where box is a 1x4 array [xmin ymin xmax ymax] specifying the box dimensions
%         x      is an Nx2 array specifying the locations of the generators
%         w      is an Nx1 array specifying the weights
%         per_x  is a boolean variable, if 'true' the computed diagram is periodic in the x-direction, if 'false' then it is not
%         per_y  is a boolean variable, if 'true' the computed diagram is periodic in the y-direction, if 'false' then it is not
%
%
%   The output is
%         area  an Nx1 array of cell areas
%         tr    an Nx1 array of cell transport costs 
%         xc    an Nx2 array of cell centroids
%         vfn   an Nx3 cell array of cell vertex coordinates, neighbours and remapped seed locations
%               provided the ith cell is non-empty then
%               vfn{i,1} contains an Nv x 2 array of vertex coordinates for the ith cell (Nv number of vertices)
%               vfn{i,2} contains an Nv x 1 array of neighbour indices for the ith cell
%               vfn{i,3} contains the remapped (into the fundamental domain) seed location
%   
%   If per_x, per_y are omitted the assumption is that the periodic flags are 'false'
%   If w is ommitted then the assumption is that the weights are all zero
%   If box is ommitted then the assumption is that box=[0 0 1 1] so the
%   domain is the unit square
%
%   It is currently possible to call mexPDall with the following arguments
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
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,per);
%     [area,tr,xc,vfn]=mexPDall_2d(x,w,per);
%
%    With four arguments, per_x=per_y=per
%
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,w,per);
%
%    With five arguments
%
%     [area,tr,xc,vfn]=mexPDall_2d(bx,x,w,per_x,per_y);
%