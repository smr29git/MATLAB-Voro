function h=plot_cells2d(vfn,varargin)

% h=plot_cells2d(vfn,OPTid,OPTcellcolour);
% 
% Plot individual grains
%   input aguments
%
%   vfn      : vertex and neighbour data
%   OPTid    : an optional list of cell ids in an nid x 1 vector
%   OPTcellcolour: an optional argument that contains either a single colour for all the cells or a
%              list of colours for each of the cells to be plotted
%    
%   output arguments
%
%   h : cell array of handles to patch objects that make up the faces of each of the plotted cells
%       this can be used to alter the styles after the plot    
%   id: the cell ids that were plotted
    
% Obtain the number of cells
    [Nc,~]=size(vfn);

    if(nargin>3)
        error('Too many input arguments\nUsage is plot_cells(vfn), plot_cells(vfn,id) or plot_cells(vfn,id,c) where c is a colour');
    elseif(nargin>1)
        % Get the cell IDs to plot
        id=varargin{1};

        % Get the number of cells to plot
        [nid,mid]=size(id);
        if(mid>1)
            error('The lits of cell ids should be an nid x 1 vector')
        end
        
        % Check that the minimum id no less than 1 and the maximum is no greater than Nc
        if(min(id)<1 || max(id)>Nc)
            error('The cell ids must be in the range [1,Nc] where Nc is the number of cells')
        end
        
        % Get unique cellids to be plotted
        [id,iid]=unique(id);

        if(nargin==3)
            % Here we assume that arguments are (vfn,id,c)
            % Get the grain colours
            cellcolour=varargin{2};
            [nc,~]=size(cellcolour);
            if(nc==1)
                cellcolour=repmat(cellcolour,nid,1);
            elseif(nc~=nid)
                error('The number of cells to plot and the length of the list of colours should be the same');
            else
                cellcolour=cellcolour(iid,:);
            end
        else
            % Here we assume that the arguments are (vfn,id)
            [nid,~]=size(id);
            cellcolour=repmat('r',nid,1);
        end
    else
        id=(1:Nc)';nid=Nc;
        cellcolour=repmat('r',nid,1);
    end

    % Get the new number of cells to plotted (after unique)
    [nid,~]=size(id);

    % Loop through each of the cells
    for j=1:nid,
        cid=id(j);
        % Can only plot non-empty cells
        if(~isempty(vfn{cid,1}))
            
            % Obtain the vertices, faces and neighbours of the cell
            verts=vfn{cid,1};
            verts=[verts;verts(1,:)];

            hold on
            h(j)=patch(verts(:,1),verts(:,2),cellcolour(j,:));
            set(h(j),'LineWidth',2);
            hold off

        end
    end

    %            plot(x(i,1),x(i,2),'r.','Markersize',5);
    
end

