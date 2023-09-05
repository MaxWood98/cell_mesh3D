%Function to import cell_mesh3d volume mesh
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.1
%Updated 01-06-2023

%Function -----------------------------------------------------------------
function [Ncell,Nface,Nvtx,faces,vertices,Nvface,cell_lr,Nface_mesh,Nface_surface] = import_mesh_cm3d(filepath)

    %Open file
    fid = fopen(filepath);
    
    %Load mesh properties 
    mdata = textscan(fid,'%d %d %d %d %d',1);
    Ncell = mdata{1,1};
    Nface = mdata{1,2};
    Nvtx = mdata{1,3};
    Nface_mesh = mdata{1,4};
    Nface_surface = mdata{1,5};

    %Load number of vertices per face
    Nvface = textscan(fid,'%d',Nface);
    Nvface = Nvface{1};
    
    %Load faces
    maxfv = max(Nvface);
    faces = zeros(Nface,maxfv);
    fbase = ' %f';
    frmat = '%f';
    for aa=1:maxfv
        frmat = strcat(frmat,fbase);
    end
    finput = textscan(fid,frmat,Nface);
    for aa=1:maxfv
        faces(:,aa) = finput{1,aa};
    end
    
    %Load adjacent cells 
    cell_lr = zeros(Nface,2);
    cell_lr_in = textscan(fid,'%d %d',Nface);
    cell_lr(:,1) = cell_lr_in{1,1};
    cell_lr(:,2) = cell_lr_in{1,2};
    
    %Load vertices
    vertices = zeros(Nvtx,3);
    vertices_in = textscan(fid,'%f %f %f',Nvtx);
    vertices(:,1) = vertices_in{1,1};
    vertices(:,2) = vertices_in{1,2};
    vertices(:,3) = vertices_in{1,3};
    
    %Close file
    fclose(fid);
end