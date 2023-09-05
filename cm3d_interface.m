%Octree Cutcell Meshing Program (Cell Mesh 3D) Interface
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.0
%Updated 24-05-2023

%Reset workspace
clearvars
clc

%TODO


%% Input

%General options 
cm3dop.condisp = 1;           %Toggle console display (1 = yes || 0 = no)
cm3dop.normdir = 1;           %Geometry surface normal direction convention (1 = normal | -1 = reverse)

%Mesh format options cm3dop
cm3dop.meshtype = 0;          %Type of mesh (0 = cutcell | 1 = minD block mesh)
cm3dop.meshinout = 'out';     %Mesh inside or outside of geometry (default out)
cm3dop.surface_dir = 'in';    %Surface normal direction switch in to / out of the mesh domain (default in)
cm3dop.boundary_dir = 'in';   %Boundary normal direction switch in to / out of the mesh domain (default in)

%Cut-Cell mesh options ====================================================
%Octree options
cm3dop.nrefine = 8;           %Maximum refinement level 
cm3dop.nrefineB = 5;          %Maximum additional refinement levels in high curvature regions
cm3dop.ncell_max = 2000000;    %Maximum number of cells
cm3dop.nrflood_i = 3;         %Refinement adjacency flooding iterations at the first refinement
cm3dop.nrflood_f = 3;         %Refinement adjacency flooding iterations at the final refinement
cm3dop.nrflood_b = 5;         %Refinement adjacency flooding iterations on boosted refinement
cm3dop.fbound = 2;            %Far field distance from object centre  
cm3dop.coffset = [0.0 1.77 0.0];   %Object/mesh centre offset (x / y / z)

%Mesh cleaning options
cm3dop.fminarea = 1e-8;        %Minimum face area as a fraction of an undeformed cell face area at each refienemnt level 
cm3dop.cminvol = 0.01;         %Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell

%Mesh geometry intersection options
cm3dop.enintmax = 10;         %Maximum number of mesh-geometry intersections on each mesh edge 
cm3dop.eintpad = 0.0;         %Edge-geometry intersection search zone padding as a fraction of the edge length 
cm3dop.int_coin_tol = 1e-8;   %Intersection co-incidence tollerance as fraction of item length  

%Surface format options
cm3dop.surftype = 0;          %Geometry surface type (0 = 'simplified' | 1 = 'exact') 
cm3dop.surfRcurvM = 2.0;      %Surface curvature multiplier
cm3dop.surfRcurvNpts = 10;    %Number of vertices used to estimate local surface curvature ==============

%Mesh smoothing options
cm3dop.nsstype = 0;           %Near surface smoothing type (0 = 'none' | 1 = 'Laplacian')
cm3dop.nsvtxflood = 2;        %Vertex selection flooding iterations from surfaces
cm3dop.nlpsmooth = 100;       %Smoothing iterations 

%ADtree options
cm3dop.adtree_spad = 0.0;     %Maximum padding size of adtree search bounding boxes as multiple of cell edge length
cm3dop.adtree_maxd = 4;       %AD tree maximum depth in terms of dimension cycles (tree is 6d)

%Gradient linking options
cm3dop.glink_con = 1;         %Construct volume to surface gradient interpolation (1 = yes | 0 = no)
cm3dop.glink_nnn = 50;        %Number of nearest neighbours to use for volume to surface gradient interpolation
cm3dop.glink_nsmooth = 4;     %Number of nearby vertices used to smooth the gradient at each surface vertex

%Boundary condition options ===============================================
%Boundary condition options 
cm3dop.set_custom_bc = 1;     %Set custom boundary conditions in specifed regions (1 = yes | 0 = no)
cm3dop.rem_ffzones = 0;       %Remove any region of the mesh connected to a far field boundary condition
cm3dop.rem_iszones = 0;       %Remove any isolated region of the mesh connected only to a wall boundary condition 

%Boundary conditions on each base mesh face
cm3dop.bc_xmin = -2;
cm3dop.bc_xmax = -2;
cm3dop.bc_ymin = -3;
cm3dop.bc_ymax = -2;
cm3dop.bc_zmin = -2;
cm3dop.bc_zmax = -2;

%Boundary condition zone bounds [xmin xmax ymin ymax zmin zmax]
BC_zones_loc = [0 0 0 0 0 0;
                0 0 0 0 0 0];

%Boundary condition zone conditions
% -3 = mass flux inflow 
% -4 = mass flux outflow
% -5 = stagnation state inflow
% -6 = back pressure outflow 
BC_zones_type = [0 ; 0];





%% Meshing function

%Write options files
write_input_file_cm3d(cm3dop);
if cm3dop.set_custom_bc == 1
    write_custom_bc_file_cm3d(BC_zones_loc,BC_zones_type);
end 

%Call meshing function 
system('cell_mesh3d');

%% Load mesh

[Ncell_m,Nface_m,Nvtx_m,faces_m,vertices_m,Nvface_m,cell_lr_m,Nface_mesh_m,Nface_surface_m] = import_mesh_cm3d('io/grid');


%% Plot mesh

%Setup figure
cla reset 
hold on

%Plot mesh (full)
% patch('vertices',vertices_m,'faces',faces_m,'EdgeAlpha',1.0,'Marker','none','facecolor',[0.8 0.8 0.8],'facealpha',0.0);
% patch('vertices',vertices_m,'faces',faces_m(33978,:),'EdgeAlpha',1.0,'Marker','*','facecolor',[0.9 0.1 0.1],'facealpha',1.0);

%Plot mesh (volume and surface)
% patch('vertices',vertices_m,'faces',faces_m(1:Nface_mesh_m,:),'EdgeAlpha',0.2,'Marker','none','facecolor',[0.8 0.8 0.8],'facealpha',0.0);
% patch('vertices',vertices_m,'faces',faces_m(Nface_mesh_m+1:Nface_m,:),'EdgeAlpha',0.5,'Marker','none','facecolor',[1.0 0.8 0.8],'facealpha',1.0);

%Read input surface file 
[~,~,vertices,connectivity] = import_cell_mesh3d_surface();

%Plot object surface 
% patch('vertices',vertices,'faces',connectivity,'FaceColor',[1.0 0.8 0.8],'EdgeColor',[1.0 0.8 0.8],'EdgeAlpha',1.0,'FaceAlpha',1.0);
% patch('vertices',vertices,'faces',connectivity,'FaceColor',[1.0 0.8 0.8],'EdgeColor',[0.2 0.2 0.2],'EdgeAlpha',1.0,'FaceAlpha',1.0);
% patch('vertices',vertices,'faces',connectivity(12049,:),'FaceColor',[1.0 0.0 0.0],'EdgeAlpha',1.0,'FaceAlpha',1.0);
% vtgt = 10724;
% plot3(vertices(vtgt,1),vertices(vtgt,2),vertices(vtgt,3),'r*')
% vtgt = 10722;
% plot3(vertices(vtgt,1),vertices(vtgt,2),vertices(vtgt,3),'g*')


%Plot boundary conditions 
% for ii=1:Nface_m
%     if cell_lr_m(ii,1) == -1 %wall
%         % patch('vertices',vertices_m,'faces',faces_m(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','c');
%     elseif cell_lr_m(ii,1) == -2 %far field
%         % patch('vertices',vertices_m,'faces',faces_m(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','b');
%     elseif cell_lr_m(ii,2) == -3 || cell_lr_m(ii,2) == -5 %inflow
%         patch('vertices',vertices_m,'faces',faces_m(ii,:),'EdgeAlpha',1.0,'Marker','none','facecolor','g');
%     elseif cell_lr_m(ii,1) == -4 || cell_lr_m(ii,1) == -6 %outflow
%         % patch('vertices',vertices_m,'faces',faces_m(ii,:),'EdgeAlpha',1.0,'Marker','none','Edgecolor','r');
%     end
% end


% %Find cell midpoints 
% fmid = zeros(1,3);
% cnf = zeros(Ncell_m,1);
% cmid = zeros(Ncell_m,3);
% for ii=1:Nface_m
%     fmid(1,1) = sum(vertices_m(faces_m(ii,1:Nvface_m(ii)),1))/double(Nvface_m(ii));
%     fmid(1,2) = sum(vertices_m(faces_m(ii,1:Nvface_m(ii)),2))./double(Nvface_m(ii));
%     fmid(1,3) = sum(vertices_m(faces_m(ii,1:Nvface_m(ii)),3))./double(Nvface_m(ii));
%     cl = cell_lr_m(ii,1);
%     cr = cell_lr_m(ii,2);
%     if cl > 0
%         cmid(cl,:) = cmid(cl,:) + fmid(:)';
%         cnf(cl) = cnf(cl) + 1;
%     end
%     if cr > 0
%         cmid(cr,:) = cmid(cr,:) + fmid(:)';
%         cnf(cr) = cnf(cr) + 1;
%     end
% end 
% cmid(:,:) = cmid(:,:)./cnf(:);
% 
% 
% %Plot cell
% ctgt = 12577;
% fmid = zeros(1,3);
% for ii=1:Nface_m
%     if cell_lr_m(ii,1) == ctgt || cell_lr_m(ii,2) == ctgt
%         % if cell_lr_m(ii,1) ~= -1 && cell_lr_m(ii,2) ~= -1
%         % if cell_lr_m(ii,1) == -1 || cell_lr_m(ii,2) == -1
% 
%             %Plot face
%             % ii
%             % cell_lr_m(ii,:)
%             patch('vertices',vertices_m,'faces',faces_m(ii,:),'FaceAlpha',0.5,'EdgeAlpha',1.0,'Marker','*','facecolor','g');
% 
%             % %Face midpoint 
%             % fmid(1) = sum(vertices_m(faces_m(ii,1:Nvface_m(ii)),1))/double(Nvface_m(ii));
%             % fmid(2) = sum(vertices_m(faces_m(ii,1:Nvface_m(ii)),2))./double(Nvface_m(ii));
%             % fmid(3) = sum(vertices_m(faces_m(ii,1:Nvface_m(ii)),3))./double(Nvface_m(ii));
%             % 
%             % %Face normal 
%             % Nf = newell_normal(Nvface_m(ii),faces_m(ii,1:Nvface_m(ii)),vertices_m);
%             % Nf = Nf./norm(Nf);
%             % 
%             % %Plot normal
%             % nl = 0.1;
%             % plot3([fmid(1) fmid(1)+nl*Nf(1)],[fmid(2) fmid(2)+nl*Nf(2)],[fmid(3) fmid(3)+nl*Nf(3)],'r');
%             % 
%             % %Adjacent cells 
%             % cl = cell_lr_m(ii,1);
%             % cr = cell_lr_m(ii,2);
%             % 
%             % %Plot cell right 
%             % if cr > 0
%             %     rvec = cmid(cr,:) - fmid(:)';
%             %     plot3([fmid(1) fmid(1)+rvec(1)],[fmid(2) fmid(2)+rvec(2)],[fmid(3) fmid(3)+rvec(3)],'b','linewidth',2);
%             % end
%             % 
%             % %Plot cell left 
%             % if cl > 0
%             %     rvec = cmid(cl,:) - fmid(:)';
%             %     plot3([fmid(1) fmid(1)+rvec(1)],[fmid(2) fmid(2)+rvec(2)],[fmid(3) fmid(3)+rvec(3)],'g','linewidth',2);
%             % end
%             % 
%             % %Face midpoint
%             % plot3(fmid(1),fmid(2),fmid(3),'r.')
% 
%             %Cell midpoint  
%             plot3(cmid(ctgt,1),cmid(ctgt,2),cmid(ctgt,3),'bo','markersize',80)
%         % end
%     end
% end 



%Debug plots ==============================================================

% %Plot face
% patch('vertices',vertices_m,'faces',faces_m(59112,:),'FaceAlpha',1.0,'EdgeAlpha',1.0,'Marker','*','facecolor','g');

%Plot vertices
% vtgt = 18678;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18681;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18682;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18679;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')

% vtgt = 18678;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18679;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18813;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')

% vtgt = 18683;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18684;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18687;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18686;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 18685;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')


% vtgt = 108074;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')
% vtgt = 108073;
% plot3(vertices_m(vtgt,1),vertices_m(vtgt,2),vertices_m(vtgt,3),'r*')


% v1 = [7.8446596085101250      -0.18750000000000006       0.80306900000000070];
% v2 = [7.8539970736976148      -0.18750000000000000       0.80292643545489684];
% v3 = [7.8685039911941024      -0.18750000000000000       0.80306900000000070];
% plot3(v1(1),v1(2),v1(3),'m.','markersize',30)
% plot3(v2(1),v2(2),v2(3),'m.','markersize',30)
% plot3(v3(1),v3(2),v3(3),'m.','markersize',30)



%Plot volume mesh surface
Nfsurf = 0;
for ff=1:Nface_m
    if cell_lr_m(ff,1) == -1 || cell_lr_m(ff,2) == -1
        Nfsurf = Nfsurf + 1;
    end
end
[~,npmaxfm] = size(faces_m);
faces_mS = zeros(Nfsurf,npmaxfm);
Nfsurf = 0;
for ff=1:Nface_m
    if cell_lr_m(ff,1) == -1 || cell_lr_m(ff,2) == -1
        Nfsurf = Nfsurf + 1;
        faces_mS(Nfsurf,:) = faces_m(ff,:);
    end
end
patch('vertices',vertices_m,'faces',faces_mS,'Marker','none','facecolor',[0.6 0.6 0.8],'facealpha',1.0,'edgecolor',[0.0 0.0 0.0],'EdgeAlpha',1.0);



% %Write surface mesh 
% vmap = zeros(Nvtx_m,1);
% Nvnew = 0;
% for ff=1:Nfsurf
%     for vv=1:npmaxfm
%         vtgt = faces_mS(ff,vv);
%         if isnan(vtgt)
%             faces_mS(ff,vv) = 0;
%             vtgt = 0;
%         end
%         if vtgt > 0 
%             if vmap(vtgt) == 0
%                 Nvnew = Nvnew + 1;
%                 vmap(vtgt) = Nvnew;
%             end
%         end
%     end
% end
% vertices_s = zeros(Nvnew,3);
% for vv=1:Nvtx_m
%     if vmap(vv) ~= 0 
%         vertices_s(vmap(vv),:) = vertices_m(vv,:);
%     end
% end
% for ff=1:Nfsurf
%     for vv=1:npmaxfm
%         vtgt = faces_mS(ff,vv);
%         if vtgt > 0 
%             faces_mS(ff,vv) = vmap(vtgt);
%         end
%     end
% end
% maxpfs = 0;
% for ff=1:Nfsurf
%     for vv=1:npmaxfm
%         vtgt = faces_mS(ff,vv);
%         if vtgt == 0 
%             if vv > maxpfs
%                 maxpfs = vv;
%             end
%             break
%         end
%     end
% end
% faces_mS = faces_mS(:,1:maxpfs);
% writematrix(faces_mS,'io/faces_surf.dat') 
% writematrix(vertices_s,'io/vertices_surf.dat') 



% vtx = load('io/vertices.dat');
% plot3(vtx(:,1),vtx(:,2),vtx(:,3),'k.')

% vtxext = load('io/vtxexternal.dat');
% plot3(vtxext(:,1),vtxext(:,2),vtxext(:,3),'ro')

% vtxp = load('io/vtx2pert.dat');
% plot3(vtxp(:,1),vtxp(:,2),vtxp(:,3),'r*')

% vtxext = load('io/ctype1.dat');
% plot3(vtxext(:,1),vtxext(:,2),vtxext(:,3),'b*')

% cotype = load('io/celloftype.dat');
% plot3(cotype(:,1),cotype(:,2),cotype(:,3),'ro')

% %Test duplicate vertices 
% hold on
% Nvtx = length(vtx);
% for ii=1:Nvtx 
%     for jj=1:Nvtx
%         if ii~= jj
%             dist = norm(vtx(ii,:) - vtx(jj,:));
%             if dist <= 0.00000001
%                 ii
%                 plot(vtx(ii,1),vtx(ii,2),'r.')
%             end
%         end 
%     end
% end

% %Test import mesh
% vertices_mt = load('io/mesh_test_vtx');
% cell_lr_m = load('io/cell_lr');
% fid = fopen('io/mesh_test_faces');
% maxfv = textscan(fid,'%d',1);
% maxfv = maxfv{1};
% fbase = ' %f';
% frmat = '%f';
% for aa=1:maxfv
%     frmat = strcat(frmat,fbase);
% end
% finput = textscan(fid,frmat);
% face_nvtx = finput{1,1};
% nface = length(face_nvtx);
% faces_mt = zeros(nface,maxfv);
% for aa=2:maxfv+1
%     faces_mt(:,aa-1) = finput{1,aa};
% end
% fclose(fid);
% Nface_m = length(faces_mt(:,1));
% 
% %Test plot mesh 
% patch('vertices',vertices_mt,'faces',faces_mt(:,:),'EdgeAlpha',0.3,'Marker','none','facecolor',[0.8 0.8 0.8],'facealpha',0.0);



% %Plot intersects
% vmsint_edge = load('io/vmsint_edge.dat');
% plot3(vmsint_edge(:,1),vmsint_edge(:,2),vmsint_edge(:,3),'b.','markersize',20)
% 
% vmsint_face = load('io/vmsint_face.dat');
% plot3(vmsint_face(:,1),vmsint_face(:,2),vmsint_face(:,3),'rp','markersize',20)
% 
% vmsint_vtx = load('io/vmsint_vtx.dat');
% plot3(vmsint_vtx(:,1),vmsint_vtx(:,2),vmsint_vtx(:,3),'go','markersize',20)




% %Find faces on vertex
% vtgtp = [0.764497 0.468618 0.109706];
% [~,maxvf] = size(faces_m);
% disp('================')
% for ii=1:Nface_m
%     for vv=1:maxvf
%         vtest = faces_m(ii,vv);
%         if vtest > 0 
%             if norm(vertices_m(vtest,:) - vtgtp) <= 0.000001
%                 ii
%                 break
%             end
%         end
%     end
% end 
% faces_m(faces_m == 0) = nan;
% patch('vertices',vertices_m,'faces',faces_m(60734,1:4),'FaceColor',[1.0 0.0 0.0],'EdgeColor',[1.0 0.0 0.0],'EdgeAlpha',1.0,'FaceAlpha',1.0,'linewidth',2);
% cell_lr_m(60734,:)
% Nvface_m(60734)


% patch('vertices',vertices_m,'faces',faces_m(60734,1:4),'FaceColor',[1.0 0.0 0.0],'EdgeColor',[1.0 0.0 0.0],'EdgeAlpha',1.0,'FaceAlpha',1.0,'linewidth',2);
% cell_lr_m(60734,:)
% Nvface_m(60734)

% patch('vertices',vertices_m,'faces',faces_m(61010,1:4),'FaceColor',[1.0 0.0 0.0],'EdgeColor',[1.0 0.0 0.0],'EdgeAlpha',1.0,'FaceAlpha',1.0,'linewidth',2);
% cell_lr_m(61010,:)
% Nvface_m(61010)



% patch('vertices',vertices_m,'faces',faces_m(138497,:),'Marker','none','facecolor',[0.6 0.6 0.8],'facealpha',1.0,'edgecolor',[0.0 0.0 0.0],'EdgeAlpha',1.0);
% a = faces_m(138497,:)


%Plot clipped face
% face_ledge = load('io/face_cvtx.dat');
% for ee=1:length(face_ledge(:,1))
%     etgt = ee;
%     plot3([face_ledge(etgt,1) face_ledge(etgt,4)],[face_ledge(etgt,2) face_ledge(etgt,5)],[face_ledge(etgt,3) face_ledge(etgt,6)],'r','linewidth',4) 
% end


% %Plot cell test
% ctgt = 59;
% for ii=1:Nface_m
%     if cell_lr_m(ii,1) == ctgt || cell_lr_m(ii,2) == ctgt
%         ii
%         disp(['cell LR = ',num2str(cell_lr_m(ii,1)),'/',num2str(cell_lr_m(ii,2))])
%         patch('vertices',vertices_m,'faces',faces_m(ii,:),'FaceAlpha',0.5,'EdgeAlpha',1.0,'Marker','*','facecolor','g');
%     end
% end 


% nl = 100;
% fmid = [0.76620153283651060       -5.1065917297084019E-002   6.8973452640973862E-002];
% Nf = [-2.8651128523818532E-004   5.8095485293972130E-004  -9.7656249999999967E-004];
% plot3([fmid(1) fmid(1)+nl*Nf(1)],[fmid(2) fmid(2)+nl*Nf(2)],[fmid(3) fmid(3)+nl*Nf(3)],'r');


%Plot cell midpoints
% cmidp = load('io/cmidp.dat');
% plot3(cmidp(:,1),cmidp(:,2),cmidp(:,3),'b*')
% ctgt = 17484;
% plot3(cmidp(ctgt,1),cmidp(ctgt,2),cmidp(ctgt,3),'r.','markersize',40)

%Plot volume mesh only edges
% voledges = load('io/voledges');
% patch('vertices',vertices_m,'faces',voledges,'FaceColor',[1.0 0.0 0.0],'EdgeColor',[1.0 0.0 0.0],'EdgeAlpha',1.0,'FaceAlpha',1.0,'linewidth',2);
 
% fedges = load('io/fedges');
% patch('vertices',vertices_m,'faces',fedges,'FaceColor',[1.0 0.0 0.0],'EdgeColor',[1.0 0.0 0.0],'EdgeAlpha',1.0,'FaceAlpha',1.0,'linewidth',2);

% f_faces = load('io/f_faces');
% patch('vertices',vertices_m,'faces',faces_m(f_faces(96),:),'Marker','none','facecolor',[0.6 0.6 0.8],'facealpha',1.0,'edgecolor',[0.0 0.0 0.0],'EdgeAlpha',1.0);

%Debug plots ==============================================================


%Format
axis square
axis equal
% axis tight
xlabel('x')
ylabel('y')
zlabel('z')
hold off

% axis([-1.0253   -0.5850   -0.2460    0.2152   -0.1717    0.2687]);
% axis([-0.7641   -0.6131   -0.0621    0.0962   -0.0040    0.1471]);
% axis([0.2038    0.3464   -0.2381   -0.1261   -0.0922    0.0198]);
% axis([0.2698    0.2983   -0.1913   -0.1689   -0.0434   -0.0211]);
% axis([0.6798    0.8454   -0.1614   -0.0308    0.0509    0.1815]);
% axis([-0.8937   -0.7458   -0.0690    0.0477   -0.0870    0.0297]);
% axis([7.7539    7.9467   -0.2717   -0.1196    0.7280    0.8800]);
% axis([4.7504    4.9850   -0.1734    0.0117   -0.8877   -0.7026]);
% axis([5.0277    5.2711    0.1853    0.3773   -0.8076   -0.6156]);


% axis([7.6817    8.0992   -0.4925   -0.0554    0.5847    1.0022]);
% axis([7.8110    7.9521   -0.2702   -0.1225    0.7376    0.8787]);
% axis([7.8423    7.8473   -0.1904   -0.1851    0.8008    0.8059]);

% axis([4.6686    4.9154   -0.0802    0.1783   -1.0464   -0.7996]);
% axis([4.7699    4.9068   -0.0097    0.1337   -0.9589   -0.8220]);

% axis([7.8419    7.8643   -0.2019   -0.1785    0.7921    0.8145]);

% axis([-0.2170   -0.0868   -0.1212    0.0152   -0.1542   -0.0240]);

% axis([7.8195    7.8819   -0.2253   -0.1599    0.7743    0.8367]);

% axis([0.1807    0.2637   -0.0075    0.0794   -0.0114    0.0716]);

% axis([0.7640    0.7647    0.4684    0.4690    0.1094    0.1100]);

% axis([0.4560    0.5155    0.1307    0.1776    0.0685    0.1155]);

%% 


function [N] = newell_normal(Nvtxf,face_vtx,vertices) 

    %Initialise 
    N= zeros(1,3);
    
    %Accumulate normal vector to the face
    for vv=1:Nvtxf
        vc = vv; 
        vn = mod(vv,Nvtxf) + 1;
        vc = face_vtx(vc);
        vn = face_vtx(vn);
        vtxC(:) = vertices(vc,:);
        vtxN(:) = vertices(vn,:);
        N(1) = N(1) - 0.5*(vtxN(3) + vtxC(3))*(vtxN(2) - vtxC(2));
        N(2) = N(2) - 0.5*(vtxN(1) + vtxC(1))*(vtxN(3) - vtxC(3));
        N(3) = N(3) - 0.5*(vtxN(2) + vtxC(2))*(vtxN(1) - vtxC(1));
    end 
end 

