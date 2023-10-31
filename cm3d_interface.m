%Octree Cutcell Meshing Program (Cell Mesh 3D) Interface
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 4.0
%Updated 31-10-2023

%Reset workspace
clearvars
clc

%TODO

%% Input

%General options 
cm3dop.condisp = 1;               %Toggle console display (1 = yes || 0 = no)
cm3dop.normdir = 1;               %Geometry surface normal direction convention (1 = normal | -1 = reverse)

%Mesh format options cm3dop
cm3dop.meshtype = 0;              %Type of mesh (0 = cutcell | 1 = minD block mesh)
cm3dop.meshinout = 'out';         %Mesh inside or outside of geometry (default out)
cm3dop.surface_dir = 'in';        %Surface normal direction switch in to / out of the mesh domain (default in)
cm3dop.boundary_dir = 'in';       %Boundary normal direction switch in to / out of the mesh domain (default in)
cm3dop.meshfrmat = 'cutcell';     %Mesh output format (cutcell / su2_dual)

%Cut-Cell mesh options ====================================================
%Octree options
cm3dop.nrefine = 10;              %Maximum refinement level 
cm3dop.nrefineB = 1;              %Maximum additional refinement levels in high curvature regions
cm3dop.ncell_max = 1000000;       %Maximum number of cells
cm3dop.nrflood_i = 2;             %Refinement adjacency flooding iterations at the first refinement
cm3dop.nrflood_f = 2;             %Refinement adjacency flooding iterations at the final refinement
cm3dop.nrflood_b = 3;             %Refinement adjacency flooding iterations on boosted refinement
cm3dop.fbound = 10.0;             %Far field distance from object centre  
cm3dop.coffset = [0.0 0.0 0.0];   %Object/mesh centre offset (x / y / z)

%Mesh cleaning options
cm3dop.fminarea = 1e-8;           %Minimum face area as a fraction of an undeformed cell face area at each refienemnt level 
cm3dop.cminvol = 0.01;            %Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell

%Mesh geometry intersection options
cm3dop.enintmax = 10;             %Maximum number of mesh-geometry intersections on each volume mesh edge 
cm3dop.int_coin_tol = 1e-12;      %Intersection co-incidence tollerance 
cm3dop.bary_loc_tol = 1e-8;       %Barycentric location tollerance 

%Surface format options
cm3dop.surftype = 0;              %Geometry surface type (0 = 'simplified' | 1 = 'exact') 
cm3dop.force_simplify = 1;        %Force simplification of all surface cells (1 = yes || 0 = no)
cm3dop.surfRcurvM = 2.0;          %Surface curvature multiplier

%Mesh smoothing options
cm3dop.nsstype = 0;               %Near surface smoothing type (0 = 'none' | 1 = 'Laplacian')
cm3dop.nsvtxflood = 2;            %Vertex selection flooding iterations from surfaces
cm3dop.nlpsmooth = 100;           %Smoothing iterations 

%ADtree options
cm3dop.adtree_spad = 0.0;         %Maximum padding size of adtree search bounding boxes as multiple of cell edge length
cm3dop.adtree_maxd = 4;           %AD tree maximum depth in terms of dimension cycles (tree is 6d)

%Gradient linking options
cm3dop.glink_con = 0;             %Construct and export volume to surface gradient interpolation (1 = yes | 0 = no)
cm3dop.glink_nnn = 20;            %Number of nearest neighbours to use for volume to surface gradient interpolation
cm3dop.glink_nsmooth = 0;         %Number of nearby vertices used to smooth the gradient at each surface vertex
cm3dop.glink_RBF_relax = 0.005;   %RBF interpolation smoothing relaxation parameter

%Boundary condition options ===============================================
%Boundary condition options 
cm3dop.set_custom_bc = 0;         %Set custom boundary conditions in specifed regions (1 = yes | 0 = no)
cm3dop.rem_ffzones = 0;           %Remove any region of the mesh connected to a far field boundary condition
cm3dop.rem_iszones = 0;           %Remove any isolated region of the mesh connected only to a wall boundary condition 

%Boundary conditions on each base mesh face
cm3dop.bc_xmin = -2;
cm3dop.bc_xmax = -2;
cm3dop.bc_ymin = -2;
cm3dop.bc_ymax = -2;
cm3dop.bc_zmin = -2;
cm3dop.bc_zmax = -2;

%Boundary condition zone bounds [xmin xmax ymin ymax zmin zmax]
BC_zones_loc = [0 0 0 0 0 0;
                0 0 0 0 0 0];

%Boundary condition zone conditions
% -1 = wall
% -2 = far field
% -3 = mass flux inflow 
% -4 = mass flux outflow
% -5 = stagnation state inflow
% -6 = freestream inflow
% -7 = back pressure outflow 
BC_zones_type = [0 ; 0];


%% Meshing function

%Write options files
write_input_file_cm3d(cm3dop);
if cm3dop.set_custom_bc == 1
    write_custom_bc_file_cm3d(BC_zones_loc,BC_zones_type);
end 

%Call meshing function 
system('cell_mesh3d mesh');
% system('cell_mesh3d project');

%% Load mesh

%Load
[Ncell_m,Nface_m,Nvtx_m,faces_m,vertices_m,Nvface_m,cell_lr_m,Nface_mesh_m,Nface_surface_m] = import_mesh_cm3d('io/grid');

%% Plot mesh

%Setup figure
cla reset 
hold on

%Plot mesh (full)
% patch('vertices',vertices_m,'faces',faces_m,'EdgeAlpha',1.0,'Marker','none','facecolor',[0.8 0.8 0.8],'facealpha',0.9);

%Read input surface file 
% [~,~,vertices,connectivity] = import_cell_mesh3d_surface();

%Plot object surface 
% patch('vertices',vertices,'faces',connectivity,'FaceColor',[0.9 0.8 1.0],'EdgeColor',[0.2 0.2 0.2],'EdgeAlpha',1.0,'FaceAlpha',1.0);

%Divide mesh to plot
Nwallface = 0;
Nintface = 0;
for ii=1:Nface_m
    if cell_lr_m(ii,1) == -1 %wall
        Nwallface = Nwallface + 1;
    else
        Nintface = Nintface + 1;
    end 
end
[~,nc] = size(faces_m);
faces_wall = zeros(Nwallface,nc);
faces_int = zeros(Nintface,nc);
Nwallface = 0;
Nintface = 0;
for ii=1:Nface_m
    if cell_lr_m(ii,1) == -1 %wall
        Nwallface = Nwallface + 1;
        faces_wall(Nwallface,:) = faces_m(ii,:);
    else
        Nintface = Nintface + 1;
        faces_int(Nintface,:) = faces_m(ii,:);
    end 
end
patch('vertices',vertices_m,'faces',faces_wall,'EdgeAlpha',1.0,'Marker','none','facecolor',[1.0 0.8 0.8],'facealpha',1.0); %Surface
% patch('vertices',vertices_m,'faces',faces_int,'EdgeAlpha',1.0,'Marker','none','facecolor',[0.8 0.8 0.8],'facealpha',0.0); %Internal

% %Plot surface and volume gradients 
% grad_surf = load('io/gradient_surf.dat');
% [~,~,vertices,connectivity] = import_cell_mesh3d_surface();
% quiver3(vertices(:,1),vertices(:,2),vertices(:,3),grad_surf(:,1),grad_surf(:,2),grad_surf(:,3),0,'r','maxheadsize',0.01)
% grad_vol = load('io/gradient.dat');
% quiver3(vertices_m(:,1),vertices_m(:,2),vertices_m(:,3),grad_vol(:,1),grad_vol(:,2),grad_vol(:,3),0,'b','maxheadsize',0.01)

%Format
axis square
axis equal
% axis tight
xlabel('x')
ylabel('y')
zlabel('z')
hold off