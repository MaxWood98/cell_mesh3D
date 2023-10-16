%Mesh Valence Calculation Function
%Max Wood 2022
%University of Bristol - Department of Aerospace Engineering

%Version 2.0
%Updated 12-05-22

%Takes --------------------------------------------------------------------
% faces (Nfcs,Npmaxf) {uint32} Mesh faces (V1 || V2 || V3 || V4 || ... || Vn)

%Returns ------------------------------------------------------------------
% valence (Nvtx,1) {uint32} Vertex valence

%Function -----------------------------------------------------------------
function [valence,MaxValence] = get_valence(faces)
    
    %Number of faces and vertices
    [Nfcs,Npmaxf] = size(faces);
    Nvtx = max(faces,[],'all');

    %Case of empty mesh
    if Nvtx == 0
        valence = [];
        MaxValence = 0;
        return
    end

    %Initialse valence array
    valence = zeros(Nvtx,1,'int32');

    %Upper bound of maximum valence 
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

            %Increment valence on each vertex
            valence(v1) = valence(v1) + 1;
            valence(v2) = valence(v2) + 1;
        end
    end
    UBmax_vlnc = max(valence(:));

    %Construct actual valence of each vertex
    vconnect = zeros(Nvtx,UBmax_vlnc,'int32');
    valence(:) = 0;
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
            for vv=1:UBmax_vlnc
                if vconnect(v1,vv) == v2
                    evalid = 0;
                    break
                end
                if vconnect(v2,vv) == v1
                    evalid = 0;
                    break
                end
            end

            %Add valence if new edge
            if evalid == 1
                
                %Increment valence on each vertex
                valence(v1) = valence(v1) + 1;
                valence(v2) = valence(v2) + 1;

                %Update vconnect
                for vv=1:UBmax_vlnc
                    if vconnect(v1,vv) == 0
                        vconnect(v1,vv) = v2;
                        break
                    end
                end
                for vv=1:UBmax_vlnc
                    if vconnect(v2,vv) == 0
                        vconnect(v2,vv) = v1;
                        break
                    end
                end
            end
        end
    end
    
    %Set maximum valence
    MaxValence = max(valence(:));
end
