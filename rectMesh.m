%% function [mesh, coord, Edof]=rectMesh(xmin, xmax, ymin, ymax, nelx, nely)
%----------------------------------------------------------------------------
% INPUT:
%          xmin,xmax: min/max x-coordinates of planar plate in xy-plane
%          ymin,ymax: min/max y-coordinates of planar plate in xy-plane
%
%          nelx: number of elements along the x-direction 
%          nely: number of elements along the y-direction 
%
% OUTPUT:          
%          
%          mesh: array where each column corresponds to an element and  
%                where the rows indicates the nodes associated with that 
%                element.
%
%          coord: Nnode x 2 array containing the x and y coordinates of
%                 each node (Nnode = total number of nodes)
%
%          Edof: Element topology matrix according to CALFEM
%                Each row correpsonds to an element and for that row
%                columns 2:end gives the associated degrees-of-freedom (5
%                dofs per node)
%
%          Edof = [ 1, u1, v1, w1, thetay1, thetax1, u2, v2, w2, thetay2,thetax2, ..., u4, v4, w4, thetay4, thetax4
%                   2, u1, v1, ....                                                              , thetax4
%                   .,  .,  ., ....                                                              ,   .
%                   .,  .,  ., ....                                                              ,   .
%                   nel, u1, v1, ....                                                         thetay4, thetax4];
%
%----------------------------------------------------------------------------

function [mesh, coord, Edof]=rectMesh(xmin, xmax, ymin, ymax, nelx, nely)
j=1;
k=1;
for i=1:nelx*nely
    mesh(:,i)=[(j-1)+nelx*j-(nelx-k);(j-1)+nelx*j+1-(nelx-k);(nelx+1)*(j+1)-(nelx-k);(nelx+1)*(j+1)-1-(nelx-k)];
    k=k+1;
    if (k>nelx)
        k=1;
        j=j+1;
        if(j>nely)
            break
        end
    end
end

[X,Y]=meshgrid(linspace(xmin,xmax,nelx+1),linspace(ymin,ymax,nely+1));
X=X';
Y=Y';
coord=[X(:),Y(:)];
nodedof=5;
Edof=zeros(nelx*nely,nodedof*4+1);
Edof=sparse(Edof);
Edof(:,1)=[[1:1:(nelx)*(nely)]'];
for i=1:nodedof
Edof(:,[i+1:nodedof:4*(nodedof)+1])=mesh'*nodedof-(nodedof-i);
end
