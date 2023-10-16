%BARYC Conversion
clearvars
clc
cla reset

vt1 = [0 0 0];
vt2 = [0 1 0];
vt3 = [1 0.75 0];
vtri = [vt1 ; vt2 ; vt3];
tri = [1 2 3];

xl = 1;
yl = 1;
vl1 = [xl yl 1];
vl2 = [xl yl -1];


%Base
O = vl1;
D = vl2 - vl1;
[Xi,Yi,Zi,~,~,~] = ray_plane_intersection(O,D,vt1,vt2,vt3);
vp = [Xi,Yi,Zi];
[u,v,w] = cart2baryc_onplane(vp,vt1,vt2,vt3);
u
v


disp('---------------')


%Test -------------
%Base triangle vectors 
A = vt2 - vt1;
B = vt3 - vt1;
C = vp - vt1;

%Matrix entries
M0 = dot(A,A);
M1 = dot(B,A);
M2 = M1;
M3 = dot(B,B);

%Determinant 
D = M0*M3 - M1*M2;

%Coordinates
ut = (dot(C,A)*M3 - M1*dot(C,B))/D
vt = (M0*dot(C,B) - dot(C,A)*M2)/D





hold on
patch('Faces',tri,'Vertices',vtri,'facealpha',0.5);
plot3([vl1(1),vl2(1)],[vl1(2),vl2(2)],[vl1(3),vl2(3)],'r')
% plot3(vp(1),vp(2),vp(3),'r.','markersize',20)
hold off
% axis([- 1 1 -1 1 -1 1])
axis vis3d


%% Test classify points 

ztol = 2e-2;

Npnts = 100;
x = linspace(-0.1,1.1,Npnts);
y = linspace(-0.1,1.1,Npnts);
[X,Y] = meshgrid(x,y);


hold on
for ii=1:Npnts
    for jj=1:Npnts
        vp = [X(ii,jj) , Y(ii,jj) , 0];
        [u,v,w] = get_barycentric_coordinates(vp,vt1,vt2,vt3);
        [vtx_loc] = vtx_bary_tri_loc(u,v,w,ztol);
        if vtx_loc(1) == 'v'
            plot3(vp(1),vp(2),vp(3),'r.','markersize',10)
        elseif vtx_loc(1) == 'e'
            plot3(vp(1),vp(2),vp(3),'b.','markersize',10)
        elseif vtx_loc(1) == 'i'
            plot3(vp(1),vp(2),vp(3),'g.','markersize',10)
        end
    end 
end 
hold off









%%



function [vtx_loc] = vtx_bary_tri_loc(u,v,w,ztol)

    %Determine location (bias to vertices when nearby for safety)
    vtx_loc = 'na';
    if ((u <= -ztol) || (v <= -ztol) || (w <= -ztol)) %outside triangle (in tollerance)
        vtx_loc = 'ot';
    elseif ((u >= ztol) && (v >= ztol) && (w >= ztol)) %inside triangle (in tollerance)
        vtx_loc = 'in';
    elseif ((u <= 2.0d0*ztol) && (v <= 2.0d0*ztol) && (w >= ztol)) %vertex 1 (in tollerance)
        vtx_loc = 'v1';
    elseif ((u >= ztol) && (v <= 2.0d0*ztol) && (w <= 2.0d0*ztol)) %vertex 2 (in tollerance)
        vtx_loc = 'v2';
    elseif ((u <= 2.0d0*ztol) && (v >= ztol) && (w <= 2.0d0*ztol)) %vertex 3 (in tollerance)
        vtx_loc = 'v3';
    elseif ((u >= ztol) && (v <= ztol) && (w >= ztol)) %edge (v1-v2) (in tollerance)
        vtx_loc = 'e1';
    elseif ((u >= ztol) && (v >= ztol) && (w <= ztol)) %edge (v2-v3) (in tollerance)
        vtx_loc = 'e2';
    elseif ((u <= ztol) && (v >= ztol) && (w >= ztol)) %edge (v3-v1) (in tollerance)
        vtx_loc = 'e3';
    else
        vtx_loc = 'uc';
    end  
end 





function [u,v,w] = get_barycentric_coordinates(vp,vt1,vt2,vt3)

    %Base triangle vectors 
    A = vt2 - vt1;
    B = vt3 - vt1;
    C = vp - vt1;
    
    %Matrix entries
    M0 = dot(A,A);
    M1 = dot(B,A);
    M2 = M1;
    M3 = dot(B,B);
    
    %Determinant 
    D = M0*M3 - M1*M2;
    
    %Coordinates
    u = (dot(C,A)*M3 - M1*dot(C,B))/D;
    v = (M0*dot(C,B) - dot(C,A)*M2)/D;
    w = 1 - u - v;
end







function [u,v,w] = cart2baryc_onplane(vp,vt1,vt2,vt3) 
    
    %Triangle normal 
    Nt = cross(vt2-vt1,vt3-vt1);
    
    %Set ray parameters 
    O(:) = vp(:);
    D(:) = Nt(:);
    
    %Find location on the triangle plane in barycentric coordinates
    [~,~,~,u,v,~] = ray_plane_intersection(O,D,vt1,vt2,vt3);
    w = 1 - u - v;
end 





function [Xi,Yi,Zi,u,v,ti] = ray_plane_intersection(O,D,V1,V2,V3)

    %Construct parameters
    E1(:) = V2(:) - V1(:);
    E2(:) = V3(:) - V1(:);
    T(:) = O(:) - V1(:);
    P(1) = D(2)*E2(3) - D(3)*E2(2);
    P(2) = -(D(1)*E2(3) - D(3)*E2(1));
    P(3) = D(1)*E2(2) - D(2)*E2(1);
    Q(1) = T(2)*E1(3) - T(3)*E1(2);
    Q(2) = -(T(1)*E1(3) - T(3)*E1(1));
    Q(3) = T(1)*E1(2) - T(2)*E1(1);
    
    %Identify intersection position on plane
    Mden = P(1)*E1(1) + P(2)*E1(2) + P(3)*E1(3);
    if Mden ~= 0 
        M = 1.0d0/(P(1)*E1(1) + P(2)*E1(2) + P(3)*E1(3));
        u = M*(P(1)*T(1) + P(2)*T(2) + P(3)*T(3));
        v = M*(Q(1)*D(1) + Q(2)*D(2) + Q(3)*D(3));
        ti = M*(Q(1)*E2(1) + Q(2)*E2(2) + Q(3)*E2(3));
        Xi = (1.0d0 - u - v)*V1(1) + u*V2(1) + v*V3(1);
        Yi = (1.0d0 - u - v)*V1(2) + u*V2(2) + v*V3(2);
        Zi = (1.0d0 - u - v)*V1(3) + u*V2(3) + v*V3(3);
    end 
end 