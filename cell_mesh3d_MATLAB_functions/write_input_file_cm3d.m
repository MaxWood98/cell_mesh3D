%Function to write cell_mesh3d input file 
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.2
%Updated 01-08-2023

%Function -----------------------------------------------------------------
function [] = write_input_file_cm3d(cm3dop)
fid = fopen('io\cell_mesh3d_options.dat','w+');
    fprintf(fid,'%s \n','#cell_mesh3d options file (version 0.0.1)');
    fprintf(fid,'%s \n',' ');
    
    fprintf(fid,'%s \n','#=== General Options =================================');
    fprintf(fid,'%s \n','#Toggle console display (1 = yes || 0 = no)');
    fprintf(fid,'%d \n',cm3dop.condisp);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Geometry surface normal direction convention (1 = normal | -1 = reverse)');
    fprintf(fid,'%d \n',cm3dop.normdir);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Mesh Format Options =============================');
    fprintf(fid,'%s \n','#Type of mesh (0 = cutcell | 1 = minD block mesh)');
    fprintf(fid,'%d \n',cm3dop.meshtype);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh inside or outside of geometry (-1 outside | 1 = inside)');
    fprintf(fid,'%s \n',cm3dop.meshinout);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Surface normal direction switch in to / out of the mesh domain (default in)');
    fprintf(fid,'%s \n',cm3dop.surface_dir);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Boundary normal direction switch in to / out of the mesh domain (default in)');
    fprintf(fid,'%s \n',cm3dop.boundary_dir);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Octree Options ==================================');
    fprintf(fid,'%s \n','#Maximum refinement level');
    fprintf(fid,'%d \n',cm3dop.nrefine);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum additional refinement levels in high curvature regions');
    fprintf(fid,'%d \n',cm3dop.nrefineB);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum number of cells');
    fprintf(fid,'%d \n',cm3dop.ncell_max);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the first refinement');
    fprintf(fid,'%d \n',cm3dop.nrflood_i);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the final refinement');
    fprintf(fid,'%d \n',cm3dop.nrflood_f);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations on boosted refinement');
    fprintf(fid,'%d \n',cm3dop.nrflood_b);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Far field distance from object centre');
    fprintf(fid,'%f \n',cm3dop.fbound);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Object/mesh centre offset (x / y / z)');
    fprintf(fid,'%f \n',cm3dop.coffset(1));
    fprintf(fid,'%f \n',cm3dop.coffset(2));
    fprintf(fid,'%f \n',cm3dop.coffset(3));
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Mesh Cleaning Options ===========================');
    fprintf(fid,'%s \n','#Minimum face area as a fraction of an undeformed cell face area at each refienemnt level');
    fprintf(fid,'%E \n',cm3dop.fminarea);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell');
    fprintf(fid,'%E \n',cm3dop.cminvol);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Mesh Geometry Intersection Options ==============');
    fprintf(fid,'%s \n','#Maximum number of mesh-geometry intersections on each mesh edge');
    fprintf(fid,'%d \n',cm3dop.enintmax);
    fprintf(fid,'%s \n',' ');     
    fprintf(fid,'%s \n','#Edge-geometry intersection search zone padding as a fraction of the edge length');
    fprintf(fid,'%E \n',cm3dop.eintpad);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Intersection co-incidence tollerance');
    fprintf(fid,'%E \n',cm3dop.int_coin_tol);
    fprintf(fid,'%s \n',' ');
        
    fprintf(fid,'%s \n','#=== Surface Format Options ==========================');
    fprintf(fid,'%s \n','#Geometry surface type (0 = simplified | 1 = exact)');
    fprintf(fid,'%d \n',cm3dop.surftype);
    fprintf(fid,'%s \n',' ');   
    fprintf(fid,'%s \n','#Surface curvature multiplier');
    fprintf(fid,'%f \n',cm3dop.surfRcurvM);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Number of vertices used to estimate local surface curvaturer');
    fprintf(fid,'%d \n',cm3dop.surfRcurvNpts);
    fprintf(fid,'%s \n',' ');  

    fprintf(fid,'%s \n','#=== Mesh Smoothing Options ==========================');
    fprintf(fid,'%s \n','#Near surface smoothing type (0 = none | 1 = Laplacian)');
    fprintf(fid,'%d \n',cm3dop.nsstype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Vertex selection flooding iterations from surfaces');
    fprintf(fid,'%d \n',cm3dop.nsvtxflood);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Smoothing iterations');
    fprintf(fid,'%d \n',cm3dop.nlpsmooth);
    fprintf(fid,'%s \n',' '); 
    
    fprintf(fid,'%s \n','#=== ADtree Options ==================================');
    fprintf(fid,'%s \n','#Padding size of adtree search bounding boxes as a multiple of cell edge length');
    fprintf(fid,'%f \n',cm3dop.adtree_spad);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#AD tree maximum depth as the number of dimension cycles (tree is 4d)');
    fprintf(fid,'%d \n',cm3dop.adtree_maxd);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Gradient Interpolation Options ==================');
    fprintf(fid,'%s \n','#Construct volume to surface gradient interpolation (1 = yes | 0 = no)');
    fprintf(fid,'%d \n',cm3dop.glink_con);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Number of nearest neighbours to use for volume to surface gradient interpolation');
    fprintf(fid,'%d \n',cm3dop.glink_nnn);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Number of nearby vertices used to smooth the gradient at each surface vertex');
    fprintf(fid,'%d \n',cm3dop.glink_nsmooth);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Boundary Condition Options ======================');
    fprintf(fid,'%s \n','#Set custom boundary conditions in specifed regions (1 = yes | 0 = no)');
    fprintf(fid,'%d \n',cm3dop.set_custom_bc);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a far field boundary condition');
    fprintf(fid,'%d \n',cm3dop.rem_ffzones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any isolated region of the mesh connected only to a wall boundary condition');
    fprintf(fid,'%d \n',cm3dop.rem_iszones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Boundary conditions on each base mesh face [xmin / xmax / ymin / ymax / zmin / zmax]');
    fprintf(fid,'%d \n',cm3dop.bc_xmin);
    fprintf(fid,'%d \n',cm3dop.bc_xmax);
    fprintf(fid,'%d \n',cm3dop.bc_ymin);
    fprintf(fid,'%d \n',cm3dop.bc_ymax);
    fprintf(fid,'%d \n',cm3dop.bc_zmin);
    fprintf(fid,'%d \n',cm3dop.bc_zmax);
    fprintf(fid,'%s \n',' '); 
fclose(fid);
end