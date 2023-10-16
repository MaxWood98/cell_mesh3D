%Mesh Connectivity Generation Function
%Max Wood 2022
%University of Bristol - Department of Aerospace Engineering

%Version 4.0
%Updated 12-05-22

%Takes --------------------------------------------------------------------
% faces (Nfcs,Npmaxf) {uint32} Mesh faces (V1 || V2 || V3 || V4 || ... || Vn)
% edges (Nedg,2) {uint32} Mesh edges (V1 || V2)
% valence (Nvtx,1) {uint32}  Mesh vertex valences
% MaxValence {double} Must be >= max(maximum valence , maximum number of points in a face) for any vertex or face in the mesh

%Returns ------------------------------------------------------------------
% V2V {uint32} Vertex to attached vertices by edges only (row is vertex number)
% V2E {uint32} Vertex to attached edges (row is vertex number)
% V2F {uint32} Vertex to attached faces (row is vertex number)
% F2E {uint32} Face to constituate edges (ordered to correct face-vertex orientation) (row is face number)
% E2F {uint32} Edge to attached faces (row is edge number)
% Npf {double} Number of vetrticies making up each mesh face (row is face number)

%Function -----------------------------------------------------------------
function [V2V,V2E,V2F,F2E,E2F,Npf] = get_connectivity(faces,edges,valence,MaxValence)

    %Number of faces edges and vertices
    [Nfcs,Npmaxf] = size(faces);
    [Nedg,~] = size(edges);
    Nvtx = max(faces,[],'all');    

    %Case of empty mesh
    if Nvtx == 0
        V2V = [];
        V2E = [];
        V2F = [];
        F2E = [];
        E2F = [];
        Npf = [];
        return
    end
    
    %Vertex to Vertex (V2V) and (Npf) =====================================
    V2V = zeros(Nvtx,MaxValence,'int32');
    Npf = zeros(Nfcs,1,'int32');
    for ff=1:Nfcs
        
        %Construct Npf in this face
        Npf_c = 0;
        for ii=1:Npmaxf
           if faces(ff,ii) ~= 0
               Npf_c = Npf_c + 1;
           end
        end
        Npf(ff) = Npf_c;
        
        %Each edge in face
        for ee=1:Npf_c

            %Edge ends 
            e1 = mod(ee-1,Npf_c) + 1;
            e2 = mod(ee,Npf_c) + 1;

            %Vertices on these ends
            v1 = faces(ff,e1);
            v2 = faces(ff,e2);

            %Check against vconnect 
            evalid = 1;
            for vv=1:MaxValence
                if V2V(v1,vv) == v2
                    evalid = 0;
                    break
                end
                if V2V(v2,vv) == v1
                    evalid = 0;
                    break
                end
            end

            %Add connection if valid (as new)
            if evalid == 1
                for vv=1:MaxValence
                    if V2V(v1,vv) == 0
                        V2V(v1,vv) = v2;
                        break
                    end
                end
                for vv=1:MaxValence
                    if V2V(v2,vv) == 0
                        V2V(v2,vv) = v1;
                        break
                    end
                end
            end
        end
    end
    
    %Vertex to Edge (V2E) =================================================
    V2E = zeros(Nvtx,MaxValence,'int32');
    for ee=1:Nedg
        
        %Vertices on this edge
        v1 = edges(ee,1);
        v2 = edges(ee,2);

        %Add to V2E
        for vv=1:MaxValence
            if V2E(v1,vv) == 0
                V2E(v1,vv) = ee;
                break
            end
        end
        for vv=1:MaxValence
            if V2E(v2,vv) == 0
                V2E(v2,vv) = ee;
                break
            end
        end
    end

    %Vertex to Face (V2F) =================================================
    V2F = zeros(Nvtx,MaxValence,'int32');
    for ff=1:Nfcs
        for ee=1:Npf(ff)
            vtxc = faces(ff,ee);
            for vv=1:MaxValence
                if V2F(vtxc,vv) == 0
                    V2F(vtxc,vv) = ff;
                    break
                end
            end
        end
    end

    %Face to Edge (F2E) ===================================================
    F2E = zeros(Nvtx,Npmaxf,'int32');
    for ff=1:Nfcs

        %Each edge in face
        f2eins = 0;
        Npf_c = Npf(ff);
        for ee=1:Npf_c

            %Edge ends 
            e1 = mod(ee-1,Npf_c) + 1;
            e2 = mod(ee,Npf_c) + 1;

            %Vertices on these ends
            v1 = faces(ff,e1);
            v2 = faces(ff,e2);

            %Find edge joining v1 v2 in this face from V2E of v1
            for vv=1:valence(v1)
                edgc = V2E(v1,vv);
                if edges(edgc,1) == v2 || edges(edgc,2) == v2
                    F2E(ff,f2eins+1) = edgc;
                    f2eins = f2eins + 1;
                    break
                end
            end
        end
    end
    
    %Edge to Face (E2F) ===================================================
    E2F = zeros(Nedg,2,'int32');
    for ff=1:Nfcs

        %Each edge in face
        Npf_c = Npf(ff);
        for ee=1:Npf_c
            edgc = F2E(ff,ee);
            if E2F(edgc,1) ~= ff && E2F(edgc,2) ~= ff
                if E2F(edgc,1) == 0
                    E2F(edgc,1) = ff;
                elseif E2F(edgc,2) == 0
                    E2F(edgc,2) = ff;
                end
            end
        end
    end
end