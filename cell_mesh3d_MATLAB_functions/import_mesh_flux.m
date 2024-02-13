%Function to import cell_mesh3d volume mesh
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.1
%Updated 01-06-2023

%Function -----------------------------------------------------------------
function [grid] = import_mesh_flux(filepath)

    %Open mesh file 
    fid = fopen(filepath,'r');
    
    %Read mesh properties
    meshdata = textscan(fid,'%d %d %d',1);
    grid.ncell = meshdata{1,1};
    grid.nface = meshdata{1,2};
    grid.nvtx = meshdata{1,3};
    
    %Build face property format
    fpr_frmt = repmat(' %d',1,grid.nface-1);
    fpr_frmt = ['%d' fpr_frmt];

    %Read face numbers of vertices 
    fnv = textscan(fid,fpr_frmt,1);
    grid.face_nvtx = zeros(grid.nface,1);
    for ii=1:grid.nface
        grid.face_nvtx(ii,1) = fnv{1,ii};
    end
    grid.face_max_nvtx = max(grid.face_nvtx(:,1));

    %Read faces 
    grid.faces = zeros(grid.nface,grid.face_max_nvtx);
    f_frmt = repmat(' %d',1,grid.face_max_nvtx-1);
    f_frmt = ['%d' f_frmt];
    face_data = textscan(fid,f_frmt,grid.nface);
    for ii=1:grid.face_max_nvtx
        grid.faces(:,ii) = face_data{1,ii};
    end

    %Read face adjacency 
    fcl = textscan(fid,fpr_frmt,1);
    fcr = textscan(fid,fpr_frmt,1);
    grid.cell_left = zeros(grid.nface,1);
    grid.cell_right = zeros(grid.nface,1);
    for ii=1:grid.nface
        grid.cell_left(ii,1) = fcl{1,ii};
        grid.cell_right(ii,1) = fcr{1,ii};
    end

    %Read vertices 
    vtxdata = textscan(fid,'%f %f %f',grid.nvtx);
    grid.vertices = zeros(grid.nvtx,3);
    for ii=1:3
        grid.vertices(:,ii) = vtxdata{1,ii};
    end

    %Close file 
    fclose(fid);
end