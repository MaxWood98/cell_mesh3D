%Function to write cell_mesh3d input file 
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 3.1
%Updated 13-02-2024

%Function -----------------------------------------------------------------
function [] = write_input_file_cm3d(cm3dopt)
fid = fopen('io\cell_mesh3d_options.dat','w+');
    fprintf(fid,'%s \n','#cell_mesh3d options file (version 0.4.2)');
    fprintf(fid,'%s \n',' ');
    
    fprintf(fid,'%s \n','#=== General Options =================================');
    fprintf(fid,'%s \n','#Toggle console display (1 = yes || 0 = no)');
    fprintf(fid,'%s %d \n','condisp =',cm3dopt.condisp);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Mesh Format Options =============================');
    fprintf(fid,'%s \n','#Type of mesh (0 = cutcell | 1 = minD block mesh)');
    fprintf(fid,'%s %d \n','meshtype =',cm3dopt.meshtype);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh inside or outside of geometry');
    fprintf(fid,'%s %s \n','meshinout =',cm3dopt.meshinout);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Surface normal direction switch in to / out of the mesh domain');
    fprintf(fid,'%s %s \n','surfnormdir =',cm3dopt.surface_dir);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Boundary normal direction switch in to / out of the mesh domain');
    fprintf(fid,'%s %s \n','bndrynormdir =',cm3dopt.boundary_dir);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Mesh output format (cutcell / su2_cutcell / su2_dual)');
    fprintf(fid,'%s %s \n','meshformat =',cm3dopt.meshfrmat);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Octree Options ================================');
    fprintf(fid,'%s \n','#Maximum refinement level');
    fprintf(fid,'%s %d \n','notrefine =',cm3dopt.nrefine);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum additional refinement levels in high curvature regions');
    fprintf(fid,'%s %d \n','nboostotrefine =',cm3dopt.nrefineB);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Maximum number of cells');
    fprintf(fid,'%s %d \n','ncellmax =',cm3dopt.ncell_max);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the first refinement');
    fprintf(fid,'%s %d \n','nadjfloodi =',cm3dopt.nrflood_i);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations at the final refinement');
    fprintf(fid,'%s %d \n','nadjfloodf =',cm3dopt.nrflood_f);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Refinement adjacency flooding iterations on boosted refinement');
    fprintf(fid,'%s %d \n','nadjfloodb =',cm3dopt.nrflood_b);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Far field distance from object centre');
    fprintf(fid,'%s %f \n','farfielddist =',cm3dopt.fbound);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Object/mesh centre offset (x / y / z)');
    fprintf(fid,'%s %f \n','offsett_x =',cm3dopt.coffset(1));
    fprintf(fid,'%s %f \n','offsett_y =',cm3dopt.coffset(2));
    fprintf(fid,'%s %f \n','offsett_z =',cm3dopt.coffset(3));
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Custom Domain Bound Specification  ==============');
    fprintf(fid,'%s \n','#Enable custom mesh domain bounds (1 = yes | 0 = no)');
    fprintf(fid,'%s %d \n','forcebounds = ',cm3dopt.set_mbounds);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Mesh domain bounds (xmin xmax ymin ymax zmin zmax)');
    fprintf(fid,'%s %E \n','bound_xmin =',cm3dopt.xmin);
    fprintf(fid,'%s %E \n','bound_xmax =',cm3dopt.xmax);
    fprintf(fid,'%s %E \n','bound_ymin =',cm3dopt.ymin);
    fprintf(fid,'%s %E \n','bound_ymax =',cm3dopt.ymax);
    fprintf(fid,'%s %E \n','bound_zmin =',cm3dopt.zmin);
    fprintf(fid,'%s %E \n','bound_zmax =',cm3dopt.zmax);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Mesh Cleaning Options ===========================');
    fprintf(fid,'%s \n','#Minimum face area as a fraction of an undeformed cell face area at each refienemnt level ');
    fprintf(fid,'%s %E \n','fminarea =',cm3dopt.fminarea);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Volume fraction of an undeformed cell at each refinement level below which a cell is classed as a sliver cell');
    fprintf(fid,'%s %E \n','cminvol =',cm3dopt.cminvol);
    fprintf(fid,'%s \n',' '); 

    fprintf(fid,'%s \n','#=== Mesh Geometry Intersection Options ==============');
    fprintf(fid,'%s \n','#Maximum number of mesh-geometry intersections on each mesh edge');
    fprintf(fid,'%s %d \n','enintmax =',cm3dopt.enintmax);
    fprintf(fid,'%s \n',' ');     
    fprintf(fid,'%s \n','#Intersection co-incidence tollerance as fraction of item length');
    fprintf(fid,'%s %E \n','intcointol =',cm3dopt.int_coin_tol);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Barycentric co-incidence tollerance on the unit triangle');
    fprintf(fid,'%s %E \n','barycointol =',cm3dopt.bary_coin_tol);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Octree cell overlap padding distance as a multiple of the cell half diagonal length');
    fprintf(fid,'%s %E \n','otrcellpad =',cm3dopt.otrcpad);
    fprintf(fid,'%s \n',' ');
        
    fprintf(fid,'%s \n','#=== Surface Format Options ==========================');
    fprintf(fid,'%s \n','#Geometry surface type (0 = simplified | 1 = exact)');
    fprintf(fid,'%s %d \n','surftype =',cm3dopt.surftype);
    fprintf(fid,'%s \n',' ');  
    fprintf(fid,'%s \n','#Force simplification of all surface cells (1 = yes || 0 = no)');
    fprintf(fid,'%s %d \n','forcesimplify =',cm3dopt.force_simplify);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Surface curvature multiplier');
    fprintf(fid,'%s %f \n','scurvmult =',cm3dopt.surfRcurvM);
    fprintf(fid,'%s \n',' ');  

    fprintf(fid,'%s \n','#=== Mesh Smoothing Options ==========================');
    fprintf(fid,'%s \n','#Near surface smoothing type (0 = none | 1 = Laplacian)');
    fprintf(fid,'%s %d \n','smoothingtype = ',cm3dopt.nsstype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Vertex selection flooding iterations from surfaces');
    fprintf(fid,'%s %d \n','smoothingnflood = ',cm3dopt.nsvtxflood);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Smoothing iterations');
    fprintf(fid,'%s %d \n','smoothingniter =',cm3dopt.nlpsmooth);
    fprintf(fid,'%s \n',' '); 
    
    fprintf(fid,'%s \n','#=== ADtree Options ==================================');
    fprintf(fid,'%s \n','#Padding size of adtree search bounding boxes as a multiple of cell edge length');
    fprintf(fid,'%s %f \n','adtpadding =',cm3dopt.adtree_spad);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#AD tree maximum depth as the number of dimension cycles (tree is 6d)');
    fprintf(fid,'%s %d \n','adtndimcycle =',cm3dopt.adtree_maxd);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#AD tree minimum divisible node size');
    fprintf(fid,'%s %d \n','adtmindivnsize =',cm3dopt.adtree_mindivnsize);
    fprintf(fid,'%s \n',' '); 
    
    fprintf(fid,'%s \n','#=== Gradient Interpolation Options ==================');
    fprintf(fid,'%s \n','#Construct volume to surface gradient interpolation (1 = yes | 0 = no)');
    fprintf(fid,'%s %d \n','glinkconstruct =',cm3dopt.glink_con);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Gradient linking type (rbf or int)');
    fprintf(fid,'%s %s \n','glinktype =',cm3dopt.glinktype);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Number of nearest neighbours to use for RBF volume to surface gradient interpolation');
    fprintf(fid,'%s %d \n','glinkrbfnpnt =',cm3dopt.glink_nnn);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#Number of vertices each side used to smooth the gradient at each surface vertex');
    fprintf(fid,'%s %d \n','glinknpntsmooth =',cm3dopt.glink_nsmooth);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Radial Basis Function Options ===================');
    fprintf(fid,'%s \n','#RBF interpolation support radius as a multiple of the maximum distance between points');
    fprintf(fid,'%s %f \n','rbfrsup =',cm3dopt.RBF_rsup);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#RBF interpolation relaxation distance as fraction of the maximum distance between points');
    fprintf(fid,'%s %f \n','rbfrrelax =',cm3dopt.RBF_relaxD);
    fprintf(fid,'%s \n',' ');
    fprintf(fid,'%s \n','#RBF interpolation relaxation parameter');
    fprintf(fid,'%s %f \n','rbfrelaxp =',cm3dopt.RBF_relaxP);
    fprintf(fid,'%s \n',' ');

    fprintf(fid,'%s \n','#=== Boundary Condition Options ======================');
    fprintf(fid,'%s \n','#Set custom boundary conditions in specifed regions (1 = yes | 0 = no)');
    fprintf(fid,'%s %d \n','setcustombcs = ',cm3dopt.set_custom_bc);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a far field boundary condition');
    fprintf(fid,'%s %d \n','remffczones =',cm3dopt.rem_ffzones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any region of the mesh connected to a non custom set boundary condition');
    fprintf(fid,'%s %d \n','remncczones =',cm3dopt.rem_nczones);
    fprintf(fid,'%s \n',' '); 
    fprintf(fid,'%s \n','#Remove any isolated region of the mesh connected only to a wall boundary condition');
    fprintf(fid,'%s %d \n','remwallonlyzones =',cm3dopt.rem_iszones);
    fprintf(fid,'%s \n',' '); 
fclose(fid);
end