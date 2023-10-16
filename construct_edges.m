%Mesh Edge Construction Function
%Max Wood 2022
%University of Bristol - Department of Aerospace Engineering

%Version 3.0
%Updated 12-05-22

%Takes --------------------------------------------------------------------
% faces (Nfcs,Npmaxf) {uint32} Mesh faces (V1 || V2 || V3 || V4 || ... || Vn)
% MaxValence {double} Must be >= max(maximum valence , maximum number of points in a face) for any vertex or face in the mesh

%Returns ------------------------------------------------------------------
% edges (Nedg,3) {uint32} Mesh edges (V1 || V2 || (Sharpness = 0))

%Function -----------------------------------------------------------------
function [edges] = construct_edges(faces,MaxValence)
    
    %Number of faces and vertices
    [Nfcs,Npmaxf] = size(faces);
    Nvtx = max(faces,[],'all');

    %Case of empty mesh
    if Nvtx == 0
        edges = [];
        return
    end

    %Construct edges
    vconnect = zeros(Nvtx,MaxValence,'int32');
    edgesC = zeros(Nvtx + Nfcs,2,'int32');
    eins = 0;
    for ff=1:Nfcs
        
        %Construct Npf in this face
        Npf = 0;
        for ii=1:Npmaxf
           if faces(ff,ii) ~= 0
               Npf = Npf + 1;
           end
        end
        
        %Each edge in face
        for ee=1:Npf

            %Edge ends 
            e1 = mod(ee-1,Npf) + 1;
            e2 = mod(ee,Npf) + 1;

            %Vertices on these ends
            v1 = faces(ff,e1);
            v2 = faces(ff,e2);

            %Check against vconnect 
            evalid = 1;
            for vv=1:MaxValence
                if vconnect(v1,vv) == v2
                    evalid = 0;
                    break
                end
                if vconnect(v2,vv) == v1
                    evalid = 0;
                    break
                end
            end

            %Add edge if valid
            if evalid == 1
                
                %Create edge
                edgesC(eins+1,1) = v1;
                edgesC(eins+1,2) = v2;

                %Update vconnect
                for vv=1:MaxValence
                    if vconnect(v1,vv) == 0
                        vconnect(v1,vv) = v2;
                        break
                    end
                end
                for vv=1:MaxValence
                    if vconnect(v2,vv) == 0
                        vconnect(v2,vv) = v1;
                        break
                    end
                end

                %Increment edge
                eins = eins + 1;
            end
        end
    end

    %Return sized edges array
    edges = zeros(eins,2,'int32');
    edges(:,:) = edgesC(1:eins,:);
end
