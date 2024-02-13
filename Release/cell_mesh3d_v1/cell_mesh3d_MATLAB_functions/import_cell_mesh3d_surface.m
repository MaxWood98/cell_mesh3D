%Function to read cell_mesh3d surface file 
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.1
%Updated 19-12-2023

%Function -----------------------------------------------------------------
function [Nvtx,Nedge,vertices,faces] = import_cell_mesh3d_surface()

    %Read surface file
    fid = fopen('io/cell_mesh3d_surface.dat');
    surf = textscan(fid,'%f %f %f');
    fclose(fid);
    
    %Extract numver of vertices and objects
    Nvtx = surf{1}(1);
    Nedge = surf{2}(1);

    %Read vertices
    vertices = zeros(Nvtx,3);
    vertices(:,1) = surf{1}(2:Nvtx+1);
    vertices(:,2) = surf{2}(2:Nvtx+1);
    vertices(:,3) = surf{3}(2:Nvtx+1);

    %Read connectivity
    faces = zeros(Nedge,3);
    faces(:,1) = surf{1}(Nvtx+2:Nvtx+Nedge+1);
    faces(:,2) = surf{2}(Nvtx+2:Nvtx+Nedge+1);
    faces(:,3) = surf{3}(Nvtx+2:Nvtx+Nedge+1);
end