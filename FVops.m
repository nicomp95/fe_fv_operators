function [ops,mesh_reduced] = FVops(mesh)
%[ops,mesh_reduced] = FVops(mesh)
%Sets up integration, averaging and differentiation operators for
%irregular finite volume mesh based on a triangulation of 2D domains
%Christian Schoof, May 2025
%
%The code assumes a triangular mesh in two dimensions. It creates finite
%volume cells with cell "centre" at the vertices of the mesh by connecting
%the mid point of each triangle edge with the centroid of the triangle.
%
%The code creates operators on a triangle-wise basis, integrating and
%averaging over the volumes and faces created by the intersection of cell
%boundaries and cell volumes with individual triangles. That intersection
%creates three cell volumes for each triangle, associated with the three
%vertices, and three faces, associated with the boundaries separating each
%pair of vertices. Each face extends from the midpoint of the triangle edge
%connecting a pair of vertices to the centroid of the triangle.
%
%If the vertices of a triangle are numbered 1, 2 and 3, then the first face
%searates vertex 1 and 2, the second face separates vertex 1 and 3, and the
%third face separates vertex 1 and 3.
%
%The code creates operators that computes a gradient on the triangle by
%interpolating linearly between values of a function at the triangle
%vertices (and formally evaluates that gradient on each face, even though
%the gradient is constant). The code also evaluates the mean value of
%linear ("Q1") functions defined through interpolation between
%the vertices, and of quadratic ("Q2") functions defined through
%interpolation between vertices and edge mid points. These means are
%defined on each face, and on the volumes associated with vertices 1,
%2 and 3.
%
%In addition, the code evaluates the size dV of each of the volumes, and
%the product of normalvector with face area (dS n_x, dS n_y) for each of
%the faces internal to a triangle, and defines averaging operators and face
%area * normal products (dS n_x_bdy, dS n_y_bdy) for triangle edges that
%are parts of the external boundary of the domain, with (n_x_bdy, n_y_bdy)
%being the usual outward-pointing normal
%
%The labelling convention employed is follows: For vectors v storing
%triangle-wise values of quantities evaluated at the triangle vertices, the
%value at the ith vertex of the jth triangle is v_{(j-1)*3+i}. For vectors
%q storing triangle-wise values of quantities evaluated at faces,
%q_{(j-1)*3 + i) is the value at the ith face in the jth triangle (i = 1
%for the face separating vertex 1 and 2 of the triangle, i = 2 for the face
%separating vertex 1 and 3, i = 3 for the face separating vertex 2 and 3).
%
%For Q1 functions, storage of nodal values in a vector is exactly as in the
%vector above. For Q2 functions, storage in a vector w puts the value of
%the function at the ith vertex in the jth triangle at w_{(j-1)*6+i}, and
%at the ith midpoint of a triangle edge (i = 1 for the edge between vertex
%1 and 2, i = 2 for edge between vertex 1 and 3, i = 3 for edge between
%vertex 2 and 3) is at w_{(j-1)*6 + 3 + i}. These conventions are identical
%to those used for Q1 and Q2 functions in TaylorHood_v3.m.
%
%For functions boundary edges, the notation follows the natural reduction
%of that used above to the case of two end points rather than three
%vertices. For functions evaluated on the two halves of each boundary edge,
%associated with the first and second end point, the value on the ith half
%of the jth edge is contained in a vector v at v_{(j-1)*2+i}. For Q1 andQ2
%functions, the notation again follows thatin TaylorHood_v3. Q1
%functions are defined as end point values in a vector v, the value at the
%ith end point of the jth boundary edge being at v_{(j-1)*2+i}. Q2
%functions are defined at those end points (i = 1,2) and at the mid point
%(i = 3), stored at v_{(j-1)*3 + i}
%
%map_cC maps values vertices at the vertices (or nodes) of the triangular
%global mesh to vertices 1, 2 and 3 of individual triangles, and allows
%summation over all vertices of individual triangles that are identified
%with the same node in the global mesh, as required by a finite volume
%discretization of a partial differential equation.
%
%
%Input arguments:
%   mesh: structure defining the mesh, with fields
%       connect:    n_element-by-n_node_per_triangle connectivity array. Each row defines a
%               triangle in the mesh, listing the node indices for the iith
%               triangle in positions connect(ii,1:3). connect may contain
%               more than three columns, for instance to account for
%               triangle midpoint nodes
%       location:   n_nodes-by-2 node location array. The ith row gives the
%               (x,y)-coordiantes of the node with node index i.
%       n_vertex:   Number of triangle vertices in the mesh (number of
%               finite volume cells)
%       n_elements: Number of triangles in the triangulation
%       dimension:  Dimensionality of the domain (2, at present --- aim to
%               expand code to cover 1 (presumably with limited use) and 3)
%       connect_bdy:    Optional n_elements_bdy-by-n_node_per_edge
%               connectivity array for boundary edge. Each row
%               corresponds to one boundary element, giving (in order) the
%               node indices of the two end points, followed by any other
%               nodes associated with the boundary
%       n_elements_bdy: If connect_bdy is supplied, n_elements_bdy is the
%               number of boundary edges
%       matchlist: optional n_match-by-2 list of pairs of periodically matched node
%               indices, in order to impose periodicity on the basis
%               functios. Any pairs containing indices greater than
%               n_vertex are ignored on the basis that they match nodes
%               that are not triangle vertices
%
%Output:
%   ops:    operator structure with fields
%       map_cC:  3*n_elements-by-n_vertex matrix. If F is an n_vertex-by-1
%               vector of values of a function at the nodes of the
%               triagulation, then f = map_cC*F produces an
%               3*n_elements-by-1 vector f of values at the vertices of
%               each triangle, labelled as described above (the vertex
%               value at the ith vertex in the jth triangle is
%               f_{(j-1)*3+i}}.Converesely, if g is a vector of values at
%               the vertices of each triangle, F = map_cC.'*f sums over the
%               values at these triangle vertices that correspond to the
%               same node in the triangulation
%       grad_x_fc: 3*n_elements-by-3*n_elements matrix; evaluates the
%               x-derivative of a function defined through values f at the
%               vertices of each triangle on each triangle-internal face as
%               grad_x_fc*f
%       grad_y_fc: As grad_x_fc, but for y-derivatives
%       mean_fQ1: 3*n_elements-by-3*n_elements matrix, evaluates the mean
%               value of a Q1 function f on the faces in each triangle as
%               mean_fQ1*f
%       mean_fQ2: 3*n_elements-by-6*n_elements matrix, same as mean_fQ1 but
%               for Q2 functions
%       mean_cQ1: 3*n_elements-by-3*n_elements matrix, evaluates the mean
%               value of a Q1 function f on the volumes in each triangle as
%               mean_fQ1*f
%       mean_cQ2: 3*n_elements-by-6*n_elements matrix, as mean_cQ1 but for
%               Q2 functions
%       nxdS_cf: 3*n_elements-by-3*n_elements matrix. For vertex i in
%               triangle j, nxdS_cf_{(j-1)*3+i,(j-1)*3+k} is the length of
%               the face separating vertex i from vertex k in the same
%               triangle, multiplied by the x-component of the normal to
%               that face, pointing away from vertex i. All entries
%               in nxdS_cf that do not correspond to a face between two
%               different vertices in the same triangle are zero
%       nydS_cf:    Same as nxdS_cf but the x-component of the normal is
%               replaced by the y-component
%       dV_cc:  Diagonal 3*n_elements-by-3*n_elements matrix, with entries
%               giving the cell volume associated with each vertex of each
%               triangle in the domain
%       map_cC_bdy: 2*n_elements_bdy-by-n_vertex matrix, maps nodal values
%               F on the triangulation of the domain onto end points of
%               each boundary element
%       interp_Q1c: interpolation from cell centre values to Q1 nodal values
%               This is a 3*n_elements-by-3*n_elements identity matrix, but
%               included for ease of changing from Q1 to Q2 basis functions
%       interp_Q2c: 6*n_elements-by-3*n_elements matrix providing
%               interpolation from cell centre values to Q2 nodal values.
%               If f is a 3*n_elements-by-1 vector of triangle vertex
%               values of some function, then f2 = interp_Q2c gives the
%               result of linear interpolation onto the nodal values of the
%               same function expressed in terms of a Q2 basis.
%       mean_fQ1_bdy: 2*n_elements_bdy-by2*n_elements_bdy matrix.
%               mean_fQ1*f evaluates the mean of a Q1 function f over each
%               half of a boundary edge, associated with the 1st and 2nd
%               end points listed in connect_bdy
%       mean_fQ2_bdy: 2*n_elements_bdy-by3*n_elements_bdy matrix, works the
%               same as mean_fQ1_bdy but for Q2 functions
%       nxdS_cf_bdy: diagonal 2*n_elements_bdy-by-2*n_elements_bdy matrix,
%               each entry giving the product of half the edge length of
%               the corresponding boundary edge endpoint with the
%               x-component of the outward-pointing unit vector
%       nydS_cf_bdy: as nxdS_cf but for the y-component of the
%               outward-pointing unit vector
%       dS_cf_bdy: as nxdS_cf, but omitting multiplication with the normal
%               vector altogether
%       map_CCreduced:  if matchlist is supplied as an input argument, the code
%               identies all sets of mutually matched nodes in the triangulation and maps them to
%               the node in that set with the lowest index. This creates a
%               reduced set of nodes, and the code then renumbers all nodes to
%               give new mode indices ranging from 1 to n_nodes-n_rep. The
%               mapping from new to old indices is map_CCreduced: A set of
%               nodal values F on the "old" indices (before mapping to the
%               reduced set) can be obtained from nodal values at the
%               newly-labelled and reduced nodes f through F =
%               map_CCreduced*f. Identical to mapP1P1reduced for the same mesh
%               as computed by TaylorHood_v3.m
%   mesh_out: postprocessed input structure 'mesh', with same fields as
%       mesh, but the following have been altered / added
%       n_vertex: number of relabelled vertices
%       n_rep: number of redundant vertices in original mesh (this is
%           n_vrep in the output of TaylorHood_v3
%       location: n_vertex-n_rep-by-dimension location array for renumbered mesh
%           vertices, with only one from each mutually matched set of
%           vertices retained
%
%Christian Schoof, May 2025
%untested as of May 26th, 2025

%extract basic mesh information
n_vertex = mesh.n_vertex;
n_elements = mesh.n_elements;   %number of elements
dimension = mesh.dimension;     %dimensionality of domain
x = mesh.location;              %array of node locations, size n_elements-by-dimension
connect = mesh.connect;         %n_lements-by-(dimension+1) array listing the nodes connect(i,1), connect(i,2) etc that are vertices of each element labelled i
if isfield(mesh,'connect_bdy')
    n_elements_bdy = mesh.n_elements_bdy;   %number of boundary faces
    connect_bdy = mesh.connect_bdy; %n_lements-by-dimension array listing the nodes connect_bdy(i,1), connect(i,2) etc that are vertices of each boundary face labell
    normal = mesh.normal;
end
if isfield(mesh,'matchlist')
    %mapping of indices in volume
    matchlist = mesh.matchlist;
    matchlist(matchlist(:,1)>n_vertex|matchlist(:,2)>n_vertex,:) = [];  %exclude any matched points that are not vertices
    volmap = 1:n_vertex;                                   %list of all indices
    volmap(matchlist(:,1)) = matchlist(:,2);            %identify indices using old labels
    volmat = sparse((1:n_vertex).', volmap,ones(size(volmap)),n_vertex,n_vertex);   %matrix version of mapping
    %identify larger sets of matched nodes (this likely does not do
    %anything for dimensions less than three)
    volmat = volmat + volmat.';
    volmat(volmat>0) = 1;
    %iterate, extending connectivity to multistep paths until all paths are
    %identified
    tic
    disp('computation of matched sets')
    volmat0 = sparse(size(volmat,1),size(volmat,2));
    while any(any(volmat0~=volmat))
        volmat0 = volmat;
        volmat = volmat*volmat;
        volmat(volmat > 0) = 1;
    end
    toc
    %identify sets of mutually matched points
    matchsets = unique(volmat,'rows');
    %associate each node in mutually matched set with first node in set
    volmap2 = 1:n_vertex;
    disp('start identification of first member of matched sets')
    tic
    matchcount = sum(volmat,2);
    matchsets = unique(volmat(matchcount>1,:),'rows');
    for ii=1:length(matchsets(:,1))
        jj = find(matchsets(ii,:),1,'first');
        volmap2(matchsets(ii,:).'==1) = jj;
    end
    toc
    %recreate matrix version of matchlist mapping, but to unique node for
    %any larger set of mutually matched nodes
    volmat2 =  sparse((1:n_vertex).', volmap2,ones(size(volmap)),n_vertex,n_vertex);
    nodes_reduced = unique(volmap2);             %nodes that are left after doing so
    n_rep = n_vertex - length(nodes_reduced);    %number of replicated nodes
    volind = sparse(nodes_reduced,(1:n_vertex-n_rep).',ones(n_vertex-n_rep,1),n_vertex,n_vertex-n_rep);    %parse redundant indices and map from old non-redundant index labels to new indicex labels
    n1reducedton1 = volmat2*volind;                          %operator transformation
    %collate output
    mesh_reduced = mesh;
    mesh_reduced.connect = mesh_reduced.connect(:,1:3);
    mesh_reduced.n_rep = n_rep;
    mesh_reduced.n_vertex = n_vertex-n_rep;
    mesh_reduced.location = x(unique(volmap),:);
    if isfield(mesh_reduced,'connect_bdy')
        mesh_reduced.connect_bdy = mesh_reduced.connect_bdy(:,1:2);
    end
else
    n1reducedton1 = speye(n_vertex);
    mesh_reduced = mesh;
    mesh_reduced.n_rep = 0;
end
    
switch dimension
    case 1
        error('not coded yet')
    case 2
        %do operators defined on entire mesh first
        %node indices, in order
        connect1 = connect(:,1);
        connect2 = connect(:,2);
        connect3 = connect(:,3);
        %corresponding vertex coordinates
        x1 = x(connect1,1); x2 = x(connect2,1); x3 = x(connect3,1);
        y1 = x(connect1,2); y2 = x(connect2,2); y3 = x(connect3,2);
        %transformation to reference triange
        Jacobian = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1);   %Jacobian of canonical transformation
        absJ = abs(Jacobian);
        chiral = sign(Jacobian);
        absJdiag = spdiags(absJ,0,n_elements,n_elements);
        chiraldiag = spdiags(chiral,0,n_elements,n_elements);
        dxdu = x2-x1; dydu = y2-y1; dxdv = x3-x1; dydv = y3-y1;
        dxdudiag = spdiags(dxdu,0,n_elements,n_elements); dydudiag = spdiags(dydu,0,n_elements,n_elements); 
        dxdvdiag = spdiags(dxdv,0,n_elements,n_elements); dydvdiag = spdiags(dydv,0,n_elements,n_elements);
        dudx = (y3-y1)./Jacobian; dudy = - (x3-x1)./Jacobian; dvdx = -(y2-y1)./Jacobian; dvdy = (x2-x1)./Jacobian;
        dudxdiag = spdiags(dudx,0,n_elements,n_elements); dudydiag = spdiags(dudy,0,n_elements,n_elements);
        dvdxdiag = spdiags(dvdx,0,n_elements,n_elements); dvdydiag = spdiags(dvdy,0,n_elements,n_elements);
        %Mapping to internal vertex labelling from global cell centre
        %indexing (and addition operator for all internal vertices
        %corresponding to global cell centres)
        n1_int = 3*n_elements;
        connectint1 = [3*(0:n_elements-1).'+1, 3*(0:n_elements-1).'+2, 3*(0:n_elements-1).'+3];
        n1ton1int = sparse(reshape(connectint1,3*n_elements,1),reshape(connect(:,1:3),3*n_elements,1),ones(3*n_elements,1),n1_int,n_vertex); %maps column vector of mesh nodal values to internal vertex values for continuous functions
        %gradient operator
        du_fc = repmat([-1, 1, 0],3,1);
        dv_fc = repmat([-1, 0, 1],3,1);
        grad_x_fc = kron(dudxdiag,du_fc) + kron(dvdxdiag,dv_fc);
        grad_y_fc = kron(dudydiag,du_fc) + kron(dvdydiag,dv_fc);
        %interpolation
        interp_Q1c = kron(speye(n_elements),speye(3));
        interp_Q2c = kron(speye(n_elements),[1, 0, 0; 0, 1, 0; 0, 0, 1; 1/2, 1/2, 0; 1/2, 0, 1/2; 0, 1/2, 1/2]);
        %averaging operator on cell faces
        weights_mean_fQ2 = [-7/108, -7/108, -5/54, 19/27, 7/27, 7/27;...
            -7/108, -5/54, -7/108, 7/27, 19/27, 7/27; ...
            -5/54, -7/108, -7/108, 7/27, 7/27, 19/27];
        mean_fQ2 = kron(speye(n_elements),weights_mean_fQ2);
        weights_mean_fQ1 = [5/12, 5/12, 1/6; 5/12, 1/6, 5/12; 1/6, 5/12, 5/12];
        mean_fQ1 = kron(speye(n_elements),weights_mean_fQ1);
        %averaging operator on cell volumes
        weights_mean_cQ2 = 6*[19/648, -19/1296, -19/1296, 47/648, 47/648, 7/324;...
            -19/1296 19/648, -19/1296, 47/648, 7/324, 47/648;...
            -19/1296, -19/1296, 19/648, 7/324, 47/648, 47/648];
        mean_cQ2 = kron(speye(n_elements),weights_mean_cQ2);
        weights_mean_cQ1 = 6*[11/108, 7/216, 7/216;...
            7/216, 11/108, 7/216;...
            7/216, 7/216, 11/108];
        mean_cQ1 = kron(speye(n_elements),weights_mean_cQ1);
        %integration over faces
        nudSref = [-1/6, -1/3, 0; 1/6, 0, -1/6; 0, 1/3, 1/6];
        nvdSref = [1/3, 1/6, 0; -1/3, 0, -1/6; 0, -1/6, 1/6];
        %nudSref = -[0, 1/6, 1/3; -1/6, 0, 1/6; -1/3, -1/6, 0];
        %nvdSref = -[0, -1/3, -1/6; 1/3, 0, 1/6; 1/6, -1/6, 0];
        nxdS_cf = kron(dydudiag*chiraldiag,nudSref )+ kron(dydvdiag*chiraldiag,nvdSref);    
        nydS_cf = - kron(dxdudiag*chiraldiag,nudSref) - kron(dxdvdiag*chiraldiag,nvdSref);
        %integration over volumes
        dV_cc = kron(absJdiag,speye(3)/6);
        if isfield(mesh,'connect_bdy')
            %integration over boundary faces
            connect_bdy1 = connect_bdy(:,1); connect_bdy2 = connect_bdy(:,2);
            xbdy1 = x(connect_bdy1,1); xbdy2 = x(connect_bdy2,1);
            ybdy1 = x(connect_bdy1,2); ybdy2 = x(connect_bdy2,2);
            absJbdy = ((xbdy2-xbdy1).^2+(ybdy2-ybdy1).^2).^(1/2);
            %internal vertices
            n1bdy_int = 2*n_elements_bdy;
            connectbdyintn1 = [2*(0:n_elements_bdy-1).'+1, 2*(0:n_elements_bdy-1).'+2];
            n1ton1bdyint = sparse(reshape(connectbdyintn1,2*n_elements_bdy,1),reshape(connect_bdy(:,1:2),2*n_elements_bdy,1),ones(n1bdy_int,1),n1bdy_int,n_vertex);
            %averaging and integration operators
            mean_cQ1_bdy = kron(speye(n_elements_bdy),[3/4, 1/4; 1/4, 3/4]);
            mean_cQ2_bdy = kron(speye(n_elements_bdy),[5/12,-1/12,2/3;-1/12,5/12,2/3]);
            %integration over face
            %nxdS_cf_bdy = sum(kron(normal(:,1).*absJbdy/2,speye(2)),2);
            %nxdS_cf_bdy = spdiags(nxdS_cf_bdy, 0, length(nxdS_cf_bdy), length(nxdS_cf_bdy));
            %nydS_cf_bdy = sum(kron(normal(:,2).*absJbdy/2,speye(2)),2);
            %nydS_cf_bdy = spdiags(nydS_cf_bdy, 0, length(nydS_cf_bdy), length(nydS_cf_bdy));
            nxdS_cf_bdy = normal(:,1) .* absJbdy / 2;
            nxdS_cf_bdy = reshape([nxdS_cf_bdy.'; nxdS_cf_bdy.'], [], 1);
            nxdS_cf_bdy = spdiags(nxdS_cf_bdy, 0, length(nxdS_cf_bdy), length(nxdS_cf_bdy));
            nydS_cf_bdy = normal(:,2) .* absJbdy / 2;
            nydS_cf_bdy = reshape([nydS_cf_bdy.'; nydS_cf_bdy.'], [], 1);
            nydS_cf_bdy = spdiags(nydS_cf_bdy, 0, length(nydS_cf_bdy), length(nydS_cf_bdy));
            dS_cf_bdy = kron(absJbdy,speye(2));
        end
    case 3
        error('not coded yet')
end
    
%collate output
ops.map_CCreduced = n1reducedton1;
ops.map_cC = n1ton1int;
ops.grad_x_fc = grad_x_fc;
ops.grad_y_fc = grad_y_fc;
ops.interp_Q1c = interp_Q1c;
ops.interp_Q2c = interp_Q2c;
ops.mean_fQ1 = mean_fQ1;
ops.mean_fQ2 = mean_fQ2;
ops.mean_cQ1 = mean_cQ1;
ops.mean_cQ2 = mean_cQ2;
ops.nxdS_cf = nxdS_cf;
ops.nydS_cf = nydS_cf;
ops.dV_cc = dV_cc;
if isfield(mesh,'connect_bdy')
    ops.map_cC_bdy = n1ton1bdyint;
    ops.mean_cQ1_bdy = mean_cQ1_bdy;
    ops.mean_cQ2_bdy = mean_cQ2_bdy;
    ops.nxdS_cf_bdy = nxdS_cf_bdy;
    ops.nydS_cf_bdy = nydS_cf_bdy;
    ops.dS_cf_bdy = dS_cf_bdy;
end

end