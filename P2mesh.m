function mesh_out = P2mesh(mesh,bdyflag)
%mesh_out = P2mesh(mesh)
%   Adds edge midpoints to connectivity matrix for a 2-D triangular mesh,
%   optionally identifies boundary faces including midpoints. Allows
%   piecewise quadratic finite element operator generation.
%
%   Input arguments:
%       mesh:   structure defining triangular mesh, with fields
%           connect:    n_elements-by-3 connectivity array for mesh. The ith row lists the
%                       node indices of the vertices of the ith triangle in the mesh
%           location:   n_nodes-by-2 array. The ith row gives the (x,y)
%                       coordinates of the mesh vertex with node index i
%           n_elements: number of triangles in the mesh
%           n_vertex:    number of vertices in the mesh
%           matchlist:  optional n_match-by-2 list of pairs of periodically
%                       matched nodes. A boundary element is only recognized as
%                       matched to another boundary elements if both end points
%                       are matched
%       bdyflag:    optional boolean, triggers identification of boundary
%                       elements
%   Output argument:
%       mesh_out:   mesh structure with augmented mesh, with fields
%           connect:    n_elements-by-6 connectivity array for the mesh.
%                   The ith row lists the node indices of the three
%                   vertices of the ith triangle, followed by the node
%                   indices of the midpoints of the edges connecting
%                   vertices 1 and 2, 1 and 3, and 2 and 3 in order.
%           location:   new n_nodes-by-2 array giving (x,y) coordinates of
%                   all nodes in the mesh. Nodes with indices larger than
%                   the number of vertices are edge midpoint nodes.
%           n_elements: Unchanged from input, the number of mesh elements
%           n_nodes:    The total number of nodes after addition of
%                   midpoint nodes
%           n_vertex:   The number of triangle vertices in the mesh, equal
%                   n_nodes in the input mesh
%           connect_bdy: if bdyflag is set, then connect_bdy is output as
%                   the n_elements_bdy-by-3 connectivity array for mesh
%                   boundary elements (identified as mesh edges that are
%                   part of the boundary of only one of the mesh
%                   triangles). The ith row gives (in order) the node
%                   indices of end point 1, end point 2 and mid point of
%                   boundary element i.
%           normal: if bdyflag is set, normal is an n_elements_bdy-by2
%                   array in which the ith row gives the outward-pointing
%                   unit normal (n_x,n_y) to the ith boundary element
%           n_elements_bdy: number of boundary elements
%           matchlist:  an updated matchlist, which includes pairings of
%                   edge midpoints for pairs of edges in which both
%                   endpoints are matched to the endpoints of the other
%                   edge
%
%Christian Schoof, May 2025
%Not tested as of May 26th, 2025

dimension = mesh.dimension;
switch dimension
    case 1
        error('1D not coded yet')
    case 2
        %extract basic mesh information
        connect = mesh.connect;
        x = mesh.location;
        n_nodes = mesh.n_vertex;        %late change in notation; would have been better to use n_nodes rather than n_vertex in the code
        n_elements = mesh.n_elements;
        connect1 = connect(:,1);
        connect2 = connect(:,2);
        connect3 = connect(:,3);
        %introduce midpoint nodes
        edges = [connect1, connect2; connect1, connect3; connect2, connect3];   %all edges; a better ordering may promote a more compact sparsity pattern?
        edges = sort(edges,2);
        [edgeunique,~,edgeind] = unique(edges,'rows');
        n_mid = size(edgeunique,1);    %number of unique edge midpoints
        x_mid = (x(edgeunique(:,1),:)+x(edgeunique(:,2),:))/2;    %midpoint coordinates
        connect_mid = [edgeind(1:n_elements), edgeind(n_elements+(1:n_elements)), edgeind(2*n_elements+(1:n_elements))];    %midpoint node index for first, second and third edge of each element (connting nodes (1,2), (1,3) and (2,3))
        %Set up for construction of operators
        connect4 = connect_mid(:,1) + n_nodes;     %append midpoint indices to list of existing vertex indices
        connect5 = connect_mid(:,2) + n_nodes;
        connect6 = connect_mid(:,3) + n_nodes;
        connect_out = [connect, connect4, connect5, connect6];
        x_out = [x; x_mid];
        n_nodes_out = n_nodes + n_mid;
        %output
        mesh_out.dimension = dimension;
        mesh_out.connect = connect_out;
        mesh_out.location = x_out;
        mesh_out.n_elements = n_elements;
        mesh_out.n_vertex = n_nodes;
        mesh_out.n_nodes = n_nodes_out;    
        if nargin > 1 && bdyflag
            %construct all faces, including mid-points, in the same ordering
            %for an (ndim-1)-dimensional simplex (vertices, then midpoints)
            connect_faces = [connect1, connect2, connect4; ...
                connect1, connect3, connect5; ...
               connect2, connect3, connect6];  %simply appends mid-point label to existing eges array above to replicate structure for a 1D connectivity array with midpoints; structure for 3D mesh would be [connect1 connect2 connect3 connect5 connect6; connect1 connect2 connect4 connect5 connect7; connect1 connect3 connect4 connect6 connect7 etc] 
            oppositevertex = [connect3; connect2; connect1];
            %list edge indices for faces, identify which faces are unique by identifying faces corresponding to unique edges (or combinations, in higher dimensions)      
            faces = [edgeind(1:n_elements); edgeind(n_elements+(1:n_elements)); edgeind(2*n_elements+(1:n_elements))];  %each row has the edgeind identifiers for the edges that are part of the boundary of the face; for 3D, this would be [edgeind(1:n_elements), edgeind(n_elements+(1:n_elements)); edgeind(1:n_elements), edgeind(2*n_elements+(1:n_elements)); ...]
            %faces = sort(faces,2);  %placeholder for 3D; redundant here as only one edge is involved in each face
            [~,connectind,faceind] = unique(faces,'rows');    %the information contained here is the same as in 'edges' but this structure should be more easily adapted to higher dimensions
            facecount = histcounts(faceind,1:max(faceind)+1);
            surfind = connectind(facecount==1);
            connect_bdy_out = connect_faces(surfind,:);
            oppositevertex = oppositevertex(surfind);
            %fix orientation: the outside of the domain should be to the right
            %of the boundary element
            vec1=x(connect_bdy_out(:,2),:)-x(connect_bdy_out(:,1),:);
            vec2=x(oppositevertex,:)-x(connect_bdy_out(:,1),:);
            %flip=find(vec1(:,1).*vec2(:,2)-vec1(:,2).*vec2(:,1)>0); %old
            flip=find(vec1(:,1).*vec2(:,2)-vec1(:,2).*vec2(:,1)<0);
            connect_bdy_out(flip,[1,2])=connect_bdy_out(flip,[2,1]);
            vec1(flip,:) = -vec1(flip,:);
            %compute normal
            normal_out = [vec1(:,2), -vec1(:,1)];
            absnormal = (vec1(:,1).^2+vec1(:,2).^2).^(1/2)*[1 1];
            normal_out = normal_out./absnormal;
            %output
            mesh_out.connect_bdy = connect_bdy_out;
            mesh_out.normal = normal_out;
            mesh_out.n_elements_bdy = size(connect_bdy_out,1);
        end
    case 3
        error('3D not coded yet')
end

if isfield(mesh,'matchlist')
    matchlist = mesh.matchlist;
    edges_complete = [connect1, connect2, connect4; ...
        connect1, connect3, connect5; ...
       connect2, connect3, connect6];   %all edges with midpoints appended; in 2D this is the same as connect_faces above
    matchflag = zeros(n_nodes,1);
    matchflag(matchlist(:,1)) = 1;
    matchflag(matchlist(:,2)) = 1;  %matchflag is equal to 1 for any vertex that appears in matchlist
    edgematchflag = matchflag(edges_complete(:,1))+matchflag(edges_complete(:,2)); %count number of matched end points for each edge
    edges_matched = edges_complete(edgematchflag==2,:);  %identify edges with two matched end points
    map = (1:n_nodes).';  %map all second matched nodes onto first matched nodes
    map(matchlist(:,2)) = matchlist(:,1);
    edges_matched(:,1:2) = map(edges_matched(:,1:2));
    edges_matched(:,1:2) = sort(edges_matched(:,1:2),2);    %order end point node labels for all edges
    [~,uniqueind,matchind] = unique(edges_matched(:,1:2),'rows');   %make all matched edges map to one version through matchind   
    edges_matched_unique = edges_matched(uniqueind,:);
    midpoint = [edges_matched(:,3),edges_matched_unique(matchind,3)];  %matched midpoint paired with original midpoint
    matchlist_midpoint = midpoint(midpoint(:,1)~=midpoint(:,2),:);
    matchlist_out = [matchlist; matchlist_midpoint];
    %output
    mesh_out.matchlist = matchlist_out;
end

end

