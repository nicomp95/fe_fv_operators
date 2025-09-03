function [ops,mesh_reduced] = TaylorHood_v3(mesh)
%[ops,, mesh_reduced] = TaylorHood_v3(mesh)
%Sets up integration, averaging and differentiation operators for
%piecewise linear finite and quadratic elements in 2 dimensions
%
%The code assumes a triangular mesh in two dimensions. It uses two sets of
%basis functions: first, piecewise linear (spanning a function space "P1")
%or quadratic (spanning "P2") over the triangulation, equal to unity at one
%node and zero at all others. For P1 functions, those nodes are the
%vertices of the triangles, for P2 functions, these are the vertices and
%midpoints of the edges of triangles in the mesh. Second, the function
%uses basis functions that are zero on all but one triangle, discontinuous across
%the boundary of that one triangle, and either constant (spanning "Q0"), elinear (spanning "Q1") or
%quadratic (spanning "Q2") on that triangle, zero at all but one of the
%nodes associated with that triangle. For Q1 functions,
%those nodes are again the vertices of the triangle, for Q2 functions,
%those nodes are again the vertices and the midpoints of the triangle. (For
%Q0 functions, any point in the triangle will suffice, and the node in
%question is therefore not specified).
%
%It is possible to represent each of the functions involved through their nodal
%values, either globally on the mesh (for P1, P2) or locally as limiting
%values at the nodes, the limit being taken from within the triangle in
%which the function is non-zero (for Q1, Q2). Doing so involves a different
%number of nodes in each case: the number of vertices (P1), the number of
%vertices plus triangle edges (P2), three times the number of triangles
%(Q1) or six times the number of triangles (Q2). Moreover, it is always possible
%to represent functions in P1, P2 as superpositions of functions in Q1, Q2,
%but not vice versa.
%
%The code assumes that functions in P1 and P2 are repsented in practice
%through vectors giving nodal values at the the relevant nodes of the mesh.
%If n_nodes is the total number of nodes in the mesh, then these nodes are
%indexed consecutively such that the first n_vertex nodes correspond to
%triangle vertices, and functions in P1 are specified purely by nodal
%values at these vertices. The remaining n_nodes-n_vertex nodes are edge
%midpoints, and functions in P2 are specified by nodal values at vertices
%*and* node midpoints.
%
%For functions specified by nodal values that are "internal" to each
%triangle in the mesh, the code assumes functions in Q1 are specified by
%vectors of length 3*n_elements, n_elements being the number of triangles
%in the mesh. The jth nodal value in the ith triangle is then given by the
%(3*(i-1)+ij)th entry into the corresponding vector of nodal values. For a
%function Q2, that nodal value is found in the (6*(i-1)+jth) entry into the
%corresponding vector. For the latter, the code assumes that the first
%three nodal values for a given triangle are vertices 1, 2 and 3 for that
%triangle, followed by nodal values for the midpoints between vertices 1
%and 2, 1 and 3 and 2 and 3 in sequence. This is the same order in which
%the corresponding nodes correspond in the connect array for the element in
%question (see input below). Q0 functions are simply specified by an
%n_elements-by-1 vector giving their values element-wise.
%
%Input arguments:
%   mesh: structure defining the mesh, with fields
%       connect:    n_element-by-6 connectivity array. Each row defines a
%               triangle in the mesh, listing (in order) the node indices
%               of vertices 1, 2, 3 and then of the midpoints of edges from
%               vertex 1 to 2, 1 to 3 and 2 to 3. The mesh node numbering must
%               be such that that the index for any vertex is smaller than
%               the index for all edge midpoints
%       location:   n_nodes-by-2 node location array. The ith row gives the
%               (x,y)-coordiantes of the node with node index i.
%       n_vertex:   Number of triangle vertices in the mesh (number of
%               degrees of freedom in the set of P1 functions)
%       n_nodes:     Total number of nodes in the mesh (number of degrees of
%               freedom in the set of P2 functions)
%       n_elements: Number of triangles in the triangulation
%       dimension:  Dimensionality of the domain (2, at present --- aim to
%               expand code to cover 1 (presumably with limited use) and 3)
%       connect_bdy:    Optional n_elements_bdy-by-3 connectivity array for
%               boundary elements. Each row corresponds to one boundary
%               element, giving (in order) the node indices of the two end
%               points, followed by the mid point.
%       n_elements_bdy: If connect_bdy is supplied, n_elements_bdy is the
%               number of boundary elements
%       matchlist: optional n_match-by-2 list of pairs of periodically matched node
%               indices, in order to impose periodicity on the basis
%               functios
%           The required form of mesh can be generated from a triangular
%           mesh wtith no mid point nodes by running P2mesh.m
%
%Output:
%   ops:    operator structure with fields
%       map_Q2P2:   A P2 function can be expressed as a superposition of Q2
%               functions. map_Q2P2 maps nodal values of a P2 function onto
%               nodal values of a Q2 function. If f is the vector of P2
%               nodal values and F the vector of Q2 nodal values, then Q2 =
%               map_Q2P2*f
%       map_Q1P1: maps nodal values of a P1 function onto nodal values of a
%               Q1 function.
%       map_Q1Q0:   interpolates a Q0 function f0 onto Q1 nodal values f1
%               through f1= map_Q1Q0*f0;
%       map_Q2Q0: as for map_Q1Q0 but maps onto Q2 functions
%       map_Q2Q1: interpolates Q1 nodal values f1 onto Q2 nodal values f2
%           through f2 = map_Q2Q1*f1
%       int_Q2Q2: F.'*int_Q2Q2*G integrates the product of two Q2
%               functions, represented by nodal value vectors F and G.
%               (Note that the product of two P2 functions with nodal value
%               vectors f and g can therefore be integrated by forming
%               f.'*map_Q2P2.'*int_Q2Q2*map_Q2P2*g)
%       int_Q2Q2_interp: F.'*int_Q2Q2_interp*G approximates the integral of
%               the product of two Q2 functions, by multiplying the nodal
%               values of F and G and interpolating that product onto a Q2
%               function that is then integrated (while the product of two
%               Q2 functions is, in general, not in Q2). Can be useful in
%               forming the Jacobian of a functional defined by integration
%               using a similar interpolation (for instance by writing J =
%               ones(1,n_elements)*int_Q0Q2*func(F) where F is a vector of
%               Q2 nodal values
%       int_Q1Q2: F.'*int_Q2Q2*G integrates the product of a Q1 function
%               with nodal values F and a Q2 function with nodal values G
%       int_Q0Q2: F.'*int_Q0Q2*G integrates the product of a Q0 function
%               (with element-wise values F) and a Q2 function with nodal
%               values G
%       int_Q1Q1: F.'*int_Q1Q1*G integrates the product of two Q1 functions
%               with nodal values F and G
%       int_Q1Q1: F.'*int_Q1Q1*G integrates the product of two Q1 functions
%               with nodal values F and G
%       int_Q1Q1_interp: As for int_Q2Q2_interp, F.'*int_Q1Q1_interp*G
%               approximates the integral of the product of two Q1 functions
%               with nodal values F and G by forming the product of those
%               nodal values and interpolating onto another Q1 function
%       int_Q0Q1: F.'*int_Q0Q1*G integrates the product of a Q0 function
%               (with element-wise values F) and a Q1 function with nodal
%               values G      
%       vol_e: vector containing volumes of elements
%       Dx_Q1P2:    Dx_Q1P2*F produces the Q1 nodal values of the
%               x-derivative of a P2 function with P2 nodal values F
%       Dy_Q1P2:    Dy_Q1P2*F produces the Q1 nodal values of the
%               y-derivative of a P2 function with P2 nodal values F
%       Dx_Q0P1:    Dx_Q0P1*F produces the Q0 element values of the
%               x-derivative of a P1 function with P1 nodal values F
%       Dy_Q0P1:    Dy_Q012*F produces the Q0 element values of the
%               y-derivative of a P1 function with P1 nodal values F
%       map_Q2P2_bdy: same as mapQ2P2 but maps P2 nodal values to Q2
%               nodal values for boundary element
%       map_Q1P1_bdy: as for map_Q2P2_bdy
%       map_Q1Q0_bdy: works as map_Q1Q0 but for nodal values on boundary
%               elements
%       map_Q2Q0_bdy: again as for map_Q1Q0_bdy but maps onto Q2 functions
%       map_Q2Q1_bdy: interpolates Q1 boundary nodal values f1 onto Q2 boundary
%               nodal values f2 through f2 = map_Q2Q1_bdy*f1
%       int_Q2Q2_bdy, int_Q2Q2_interp_bdy, int_Q1Q2_bdy, int_Q0Q2_bdy,
%       int_Q1Q1_bdy, int_Q1Q1_bdy, int_Q1Q1_interp_bdy, int_Q0Q1_bdy:
%               equivalent of operators missing the "_bdy" suffix, but for
%               integration over the domain boundary
%       area_e: vector containing volumes of boundary elements
%       map_P2P2reduced:  if matchlist is supplied as an input argument, the code
%           identies all sets of mutually matched nodes and maps them to
%           the node in that set with the lowest index. This creates a
%           reduced set of nodes, and the code then renumbers all nodes to
%           give new mode indices ranging from 1 to n_nodes-n_rep. The
%           mapping from new to old indices is map_P2P2reduced: A set of
%           nodal values F on the "old" indices (before mapping to the
%           reduced set) can be obtained from nodal values at the
%           newly-labelled and reduced nodes f through F =
%           map_P2P2reduced*f
%       map_P1P1reduced: same as map_P2P2reduced, but done only for triangle
%           vertices rather than all nodes in the domain
%   mesh_out: postprocessed input structure 'mesh', with same fields as
%       mesh, but the following have been altered / added
%       n_nodes: number of relablled nodes (n_nodes-n_rep in terms of the original value n_nodes, and the number n_rep of
%               periodically matched node pairs)
%       n_vertex: number of relabelled vertices
%       n_rep:  number of redundnant nodes in original mesh
%       n_vrep: number of redundant vertices in original mesh
%       location: n_nodes-n_rep-by-dimension location array for renumbered mesh
%           vertices, with only one of each periodically matched set of nodes
%           retained (specifically,that with the lowest original node
%           index)
%
%Christian Schoof, May 2025; not tested as of May 26th, 2025

%extract basic mesh information
n_nodes = mesh.n_nodes;           %number of nodes
n_vertex = mesh.n_vertex;
n_elements = mesh.n_elements;   %number of elements
dimension = mesh.dimension;     %dimensionality of domain
x = mesh.location;              %array of node locations, size n_elements-by-dimension
connect = mesh.connect;         %n_lements-by-(dimension+1) array listing the nodes connect(i,1), connect(i,2) etc that are vertices of each element labelled i
if isfield(mesh,'connect_bdy')
    n_elements_bdy = mesh.n_elements_bdy;   %number of boundary faces
    connect_bdy = mesh.connect_bdy; %n_lements-by-dimension array listing the nodes connect_bdy(i,1), connect(i,2) etc that are vertices of each boundary face labell
end
if isfield(mesh,'matchlist')
    %mapping of indices in volume
    matchlist = mesh.matchlist;
    volmap = 1:n_nodes;                                   %list of all indices
    volmap(matchlist(:,1)) = matchlist(:,2);            %identify indices using old labels
    volmat = sparse((1:n_nodes).', volmap,ones(size(volmap)),n_nodes,n_nodes);   %matrix version of mapping
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
    disp('start identification of first member of matched sets')
    tic
    matchcount = sum(volmat,2);
    matchsets = unique(volmat(matchcount>1,:),'rows');
    size(matchsets)
    %associate each node in mutually matched set with first node in set
    volmap2 = 1:n_nodes;
    for ii=1:length(matchsets(:,1))
        jj = find(matchsets(ii,:),1,'first');
        volmap2(matchsets(ii,:).'==1) = jj;
    end
    toc
    %recreate matrix version of matchlist mapping, but to unique node for
    %any larger set of mutually matched nodes
    volmat2 =  sparse((1:n_nodes).', volmap2,ones(size(volmap)),n_nodes,n_nodes);
    nodes_reduced = unique(volmap2);             %nodes that are left after doing so
    n_rep = n_nodes - length(nodes_reduced);    %number of replicated nodes
    volind = sparse(nodes_reduced,(1:n_nodes-n_rep).',ones(n_nodes-n_rep,1),n_nodes,n_nodes-n_rep);    %parse redundant indices and map from old non-redundant index labels to new indicex labels
    n2reducedton2 = volmat2*volind;                          %operator transformation
    %repeat just for vertices
    volmap3 = volmap2(1:n_vertex);
    volmat3 = sparse((1:n_vertex).',volmap3,ones(size(volmap3)),n_vertex,n_vertex);
    nodes_reduced2 = unique(volmap3);
    n_rep2 = n_vertex - length(nodes_reduced2);
    volind2 = sparse(nodes_reduced2,(1:n_vertex-n_rep2).',ones(n_vertex-n_rep2,1),n_vertex,n_vertex-n_rep2);
    n1reducedton1 = volmat3*volind2;
    mesh_reduced = mesh;
    mesh_reduced.n_rep = n_rep;
    mesh_reduced.n_vrep = n_rep2;
    mesh_reduced.n_nodes = n_nodes-n_rep;
    mesh_reduced.n_vertex = n_vertex-n_rep2;
    mesh_reduced.location = x(unique(volmap));
else
    mesh_reduced = mesh;
    n1reducedton1 = speye(n_vertex);
    n2reducedton2 = speye(n_nodes);
    mesh_reduced.n_rep = 0;
    mesh_reduced.n_vrep = 0;
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
        %Jacobian determinant
        Jacobian = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1);   %Jacobian of canonical transformation
        absJ = abs(Jacobian);
        %Element-wise internal vertices at mesh nodes
        n1_int = 3*n_elements;
        connectint1 = [3*(0:n_elements-1).'+1, 3*(0:n_elements-1).'+2, 3*(0:n_elements-1).'+3];
        n1ton1int = sparse(reshape(connectint1,3*n_elements,1),reshape(connect(:,1:3),3*n_elements,1),ones(3*n_elements,1),n1_int,n_vertex); %maps column vector of mesh nodal values to internal vertex values for continuous functions
        n2_int = 6*n_elements;
        connectint2 = [6*(0:n_elements-1).'+1, 6*(0:n_elements-1).'+2, 6*(0:n_elements-1).'+3, 6*(0:n_elements-1).'+4, 6*(0:n_elements-1).'+5, 6*(0:n_elements-1).'+6];
        n2ton2int = sparse(reshape(connectint2,6*n_elements,1),reshape(connect,6*n_elements,1),ones(6*n_elements,1),n2_int,n_nodes); %maps column vector of mesh nodal values to internal vertex values for continuous functions
        %map from n1 to n2 or e to n1, n2
        n0ton1 = kron(speye(n_elements),ones(3,1));
        n0ton2 = kron(speye(n_elements),ones(6,1));
        n1ton2 = kron(speye(n_elements),[1, 0, 0; 0, 1, 0; 0, 0, 1; 1/2, 1/2, 0; 1/2, 0, 1/2; 0, 1/2, 1/2]);
        %weight function for integration on each element
        absJdiag = spdiags(absJ,0,n_elements,n_elements);
        %P2-P2 product nodal integration
        weights_n2n2 = [1/60, -1/360, -1/360, 0 , 0, -1/90; ...
            -1/360, 1/60, -1/360, 0, -1/90, 0;...
            -1/360, -1/360, 1/60, -1/90, 0, 0;...
            0, 0, -1/90, 4/45, 2/45, 2/45;...
            0, -1/90, 0, 2/45, 4/45, 2/45;...
            -1/90, 0, 0, 2/45, 2/45, 4/45];
        int_n2n2 = kron(absJdiag,weights_n2n2);
        %P2-P2 nodal integration by quadratic interpolation
        weights_n2n2_interp = sparse(4:6,4:6,1/6*ones(3,1),6,6);
        int_n2n2_interp = kron(absJdiag,weights_n2n2_interp);
        %P2 nodal integration
        weights_en2 = [0 0 0 1/6 1/6 1/6];
        int_en2 = kron(absJdiag,weights_en2);
        %P1-P2 nodal integration
        weights_n1n2 = [1/60, -1/120, -1/120, 1/15, 1/15, 1/30;...
            -1/120, 1/60, -1/120, 1/15, 1/30, 1/15;...
            -1/120, -1/120, 1/60, 1/30, 1/15, 1/15];
        int_n1n2 = kron(absJdiag,weights_n1n2);
        %P1-P1 exact nodal integration
        weights_n1n1 = [1/12, 1/24, 1/24; 1/24, 1/12, 1/24; 1/24, 1/24, 1/12];
        int_n1n1=  kron(absJdiag,weights_n1n1);
        %P1-P1 linear interpolation
        weights_n1n1_interp = eye(3)/6;
        int_n1n1_interp =  kron(absJdiag,weights_n1n1_interp);
        %P1 nodal integration
        weights_en1 = [1/6, 1/6, 1/6];
        int_en1 = kron(absJdiag,weights_en1);
        %element size
        vol_e = 1/2*absJ;
        %Differential operators
        %Jacobian matrix of element-wise canonical transformation
        dudx = (y3-y1)./Jacobian; dudy = - (x3-x1)./Jacobian; dvdx = -(y2-y1)./Jacobian; dvdy = (x2-x1)./Jacobian;
        dudxdiag = spdiags(dudx,0,n_elements,n_elements); dudydiag = spdiags(dudy,0,n_elements,n_elements);
        dvdxdiag = spdiags(dvdx,0,n_elements,n_elements); dvdydiag = spdiags(dvdy,0,n_elements,n_elements);
        du_n1n2 = [-3, -1, 0, 4, 0, 0;...
            1, 3, 0, -4, 0, 0;...
            1, -1, 0, 0, -4, 4];
        dv_n1n2 = [-3, 0, -1, 0, 4, 0;...
            1, 0, -1, -4, 0, 4;...
            1, 0, 3, 0, -4, 0];
        Dx_n1n2 = kron(dudxdiag,du_n1n2) + kron(dvdxdiag,dv_n1n2);
        Dy_n1n2 = kron(dudydiag,du_n1n2) + kron(dvdydiag,dv_n1n2);
        du_en1 = [-1, 1, 0];
        dv_en1 = [-1, 0, 1];
        Dx_en1 = kron(dudxdiag,du_en1) + kron(dvdxdiag,dv_en1);
        Dy_en1 = kron(dudydiag,du_en1) + kron(dvdydiag,dv_en1);
        
        %repeat for P1 elements
        %boundary operators
        if isfield(mesh,'connect_bdy')
            connect_bdy1 = connect_bdy(:,1); connect_bdy2 = connect_bdy(:,2); connect_bdy3 = connect_bdy(:,3);
            xbdy1 = x(connect_bdy1,1); xbdy2 = x(connect_bdy2,1);
            ybdy1 = x(connect_bdy1,2); ybdy2 = x(connect_bdy2,2);
            absJbdy = ((xbdy2-xbdy1).^2+(ybdy2-ybdy1).^2).^(1/2);
            %internal vertices
            n1bdy_int = 2*n_elements_bdy;
            connectbdyintn1 = [2*(0:n_elements_bdy-1).'+1, 2*(0:n_elements_bdy-1).'+2];
            n1ton1bdyint = sparse(reshape(connectbdyintn1,2*n_elements_bdy,1),reshape(connect_bdy(:,1:2),2*n_elements_bdy,1),ones(n1bdy_int,1),n1bdy_int,n_vertex);
            n2bdy_int = 3*n_elements_bdy;
            connectbdyintn2 = [3*(0:n_elements_bdy-1).'+1, 3*(0:n_elements_bdy-1).'+2, 3*(0:n_elements_bdy-1).'+3];
            n2ton2bdyint = sparse(reshape(connectbdyintn2,3*n_elements_bdy,1),reshape(connect_bdy,3*n_elements_bdy,1),ones(n2bdy_int,1),n2bdy_int,n_nodes);
            %map between e, n1, n2
            n0ton1bdy = kron(speye(n_elements_bdy),ones(2,1));
            n0ton2bdy = kron(speye(n_elements_bdy),ones(3,1));
            n1ton2bdy = kron(speye(n_elements_bdy),[1 0; 0 1; 1/2 1/2]);
            %weight function for each boundary element
            absJdiag_bdy = spdiags(absJbdy,0,n_elements_bdy,n_elements_bdy);
            %P2-P2 product nodal integration
            weights_n2n2_bdy = [2/15, -1/30, 1/15; -1/30, 215, 1/15; 1/15, 1/15, 8/15];
            int_n2n2_bdy = kron(absJdiag_bdy,weights_n2n2_bdy);
            %P2-P2 nodal integration by quadratic interpolation
            weights_n2n2_bdy_interp = spdiags([1/6, 1/6, 2/3].',0,3,3);
            int_n2n2_bdy_interp = kron(absJdiag_bdy,weights_n2n2_bdy_interp);
            %P2 nodal integration
            weights_en2_bdy = [1/6 1/6 2/3];
            int_en2_bdy = kron(absJdiag_bdy,weights_en2_bdy);
            %P1-P2 nodal integration
            weights_n1n2_bdy = [1/6, 0, 1/3; 0, 1/6, 1/3];
            int_n1n2_bdy = kron(absJdiag_bdy,weights_n1n2_bdy);
            %P1-P1 exact nodal integration
            weights_n1n1_bdy = [1/3, 1/6; 1/6, 1/3];
            int_n1n1_bdy =  kron(absJdiag_bdy,weights_n1n1_bdy);
            %P1-P1 linear interpolation
            weights_n1n1_bdy_interp = eye(2)/2;
            int_n1n1_bdy_interp =  kron(absJdiag_bdy,weights_n1n1_bdy_interp);
            %P1 nodal integration
            weights_en1_bdy = [1/2, 1/2];
            int_en1_bdy = kron(absJdiag_bdy,weights_en1_bdy);
            %surface element size
            area_e_bdy = absJbdy;
        end
   case 3
        error('not coded yet')
end
    
    %collate output
    ops.map_P2P2reduced = n2reducedton2;
    ops.map_P1P1reduced = n1reducedton1;
    ops.map_Q2P2 = n2ton2int;
    ops.map_Q1P1 = n1ton1int;
    ops.map_Q1Q0 = n0ton1;
    ops.map_Q2Q0 = n0ton2;
    ops.map_Q2Q1 = n1ton2;
    ops.int_Q2Q2 = int_n2n2;
    ops.int_Q2Q2_interp = int_n2n2_interp;
    ops.int_Q1Q2 = int_n1n2;
    ops.int_Q0Q2 = int_en2;
    ops.int_Q1Q1 = int_n1n1;
    ops.int_Q1Q1_interp = int_n1n1_interp;
    ops.int_Q0Q1 = int_en1;
    ops.vol_e = vol_e;
    ops.Dx_Q1P2 = Dx_n1n2*n2ton2int;
    ops.Dy_Q1P2 = Dy_n1n2*n2ton2int;
    ops.Dx_Q0P1 = Dx_en1*n1ton1int;
    ops.Dy_Q0P1 = Dy_en1*n1ton1int;
    if isfield(mesh,'connect_bdy')
        ops.map_Q1P1_bdy = n1ton1bdyint;
        ops.map_Q2P2_bdy = n2ton2bdyint;
        ops.map_Q1Q0_bdy = n0ton1bdy;
        ops.map_Q2Q0_bdy = n0ton2bdy;
        ops.map_Q2Q1_bdy = n1ton2bdy;
        ops.int_Q2Q2_bdy = int_n2n2_bdy;
        ops.int_Q2Q2_interp_bdy = int_n2n2_bdy_interp;
        ops.int_Q1Q2_bdy = int_n1n2_bdy;
        ops.int_Q0Q2_bdy = int_en2_bdy;
        ops.int_Q1Q1_bdy = int_n1n1_bdy;
        ops.int_Q1Q1_interp_bdy = int_n1n1_bdy_interp;
        ops.int_Q0Q1_bdy = int_en1_bdy;
        ops.area_e.bdy = area_e_bdy;
    end

end