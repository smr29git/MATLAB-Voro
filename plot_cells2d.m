function h=plot_cells2d(bx,x,vfn)
%
% h=plot_cells2d(bx,x,vfn)
%
% function to plot the cells
%
% bx is the bounding box in the form of limits [xmin ymin xmax ymax]
% x is an Nx2 array of seed locations
% w is an Nx1 array of weights
% vfn is an {N,2} cell array containing the vertices of each cell and the neighbour indices

%    xr=cell2mat(vfn(:,3));
    
    lightblue=[0.7,0.8,1.0];
    
    figure;
    clf;
    
    [N,~]=size(x);
    [Nv,~]=size(vfn);
    
    if(Nv~=N)
        error('x should be Nx2 array and vfn an {N,2} cell array');
    end
    
    hold on
    for i=1:N,
        V=vfn{i,1};
        if(~isempty(V))    
    
            V=[V;V(1,:)];
            h(i)=patch(V(:,1),V(:,2),lightblue);
            set(h(i),'LineWidth',2);
            plot(x(i,1),x(i,2),'r.','Markersize',5);
        end            
    end

    % Plot the bounding box
    bbx=[bx(1) bx(2);
         bx(1) bx(4);
         bx(3) bx(4);
         bx(3) bx(2);
         bx(1) bx(2);
         ];

    plot(bbx(:,1),bbx(:,2),'g-','LineWidth',2);
    %axis([min(bbx(:,1)) max(bbx(:,1)) min(bbx(:,2)) max(bbx(:,2))]);
    axis equal
    hold off
end

