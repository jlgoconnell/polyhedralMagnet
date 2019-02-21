
% Function to approximate force between two polyhedral magnets
% James O'Connell 13/7/18

function [F,T] = polyhedronForceOld(points1,points2,Mag1,Mag2,leverpoint)

tic;

points1 = uniquetol(points1,'ByRows',true);
points2 = uniquetol(points2,'ByRows',true);
Mag1 = Mag1/(4*pi*10^(-7));
Mag2 = Mag2/(4*pi*10^(-7));



% Initialisation:
polygons1 = {};
M1 = {};

%myFz = zeros(length(range),2);

% Define the first magnet
S1 = alphaShape(points1,Inf);
[mypolygons1,myM1] = addPolyhedron(points1,Mag1);
polygons1 = [polygons1,mypolygons1];
M1 = [M1,myM1];

% Now define the second magnet
S2 = alphaShape(points2,Inf);
M2 = Mag2;

% Set up a mesh for magnet2
S = alphaShape(points2,10);
[bf,P] = boundaryFacets(S);
realfaces = bf;
facearea = []; % The area of each triangle
centroid = mean(P);
n = []; % The normal vector of each triangle
tricentres = []; % The centre of each triangle
tricoords = []; % The vertices of each triangle:
                % [x1,x2,x3,y1,y2,y3,z1,z2,z3]

% Set up a while loop to refine mesh until convergence is achieved
Fx = Inf;
Fy = Inf;
Fz = Inf;
err = Inf;
pass = 0;
mintriarea = Inf;
maxtriarea = surfaceArea(S);
%while(err > convergencecriteria)
    pass = pass + 1;
    %meshsize = mintriarea*0.9;
    %meshsize = maxtriarea*0.01;
    meshsize = 0.00578703703703704*surfaceArea(S)
    %meshsize = 5e-6
    mintriarea = Inf;
    toobig = 1;
    bf = realfaces;
    realfaces = [];
    tricoords = [];
    tricentres = [];
    n = [];
    facearea = [];
    
    fprintf('Starting meshing...\n');

    while(toobig)
        mypts = [P(bf(1,:),1),P(bf(1,:),2),P(bf(1,:),3)];
        v1 = mypts(2,:)-mypts(1,:);
        v2 = mypts(3,:)-mypts(1,:);
        v3 = mypts(3,:)-mypts(2,:);
        longestside = max([norm(v1),norm(v2),norm(v3)]);
        myn = cross(v1,v2);
        area = 0.5*norm(myn);
        if area < mintriarea
            mintriarea = area;
        end
        if area > maxtriarea
            maxtriarea = area;
        end
        if area > meshsize% || longestside > 2*sqrt(meshsize)
            newpts = [mean(mypts(1:2,:)); mean(mypts(2:3,:)); mean(mypts([1,3],:))];
            P = [P; newpts];
            bf(end+1:end+4,:) = [bf(1,1),length(P)-2,length(P);bf(1,2),length(P)-2,length(P)-1;bf(1,3),length(P),length(P)-1;length(P),length(P)-1,length(P)-2];
            bf(1,:) = [];
        else
            realfaces = [realfaces; bf(1,:)];
            tricoords = [tricoords; mypts(:)'];
            tricentres = [tricentres; mean(mypts)];
            facearea = [facearea; area];
            bf(1,:) = [];
            myn = myn/norm(myn);
            if dot(myn,(tricentres(end,:)-centroid)) < 0
                myn = -myn;
            end
            n = [n;myn];
        end
        toobig = size(bf,1);
    end
    
    fprintf('Completed meshing. Starting force evaluation...\n');

    % Now we have n (normal vector) and tricentres (centre) for each triangle

    % Calculate the field due to magnet1 at the centre of each triangle of
    % magnet2
    MdotN = n*M2';
    x = tricentres(:,1);
    y = tricentres(:,2);
    z = tricentres(:,3);
    Bx = zeros(size(x));
    By = zeros(size(y));
    Bz = zeros(size(z));
    size(x)
    max(facearea)

    for p = 1:length(polygons1)
        for i = 1:length(x)
            if MdotN(i) ~= 0
                [Ba,Bb,Bc] = polygon_field(polygons1{p},[x(i),y(i),z(i)]',M1{p});
                Bx(i) = Bx(i) + Ba;
                By(i) = By(i) + Bb;
                Bz(i) = Bz(i) + Bc;
            end
        end
        fprintf('Completed %f%%...\n',p/length(polygons1)*100);
    end

    % Now numerically integrate for the force between the two magnets
    newx = Fx;
    newy = Fy;
    newz = Fz;
    Fx = sum(MdotN.*Bx.*facearea);
    Fy = sum(MdotN.*By.*facearea);
    Fz = sum(MdotN.*Bz.*facearea);
    B = [Bx,By,Bz];
    myleverpoint = [x,y,z]-leverpoint;
    tau = sum(MdotN.*cross(myleverpoint,B).*facearea);
    taux = tau(1);
    tauy = tau(2);
    tauz = tau(3);
    
    max(facearea)
    
    
    err = abs(max([abs(Fx-newx),abs(Fy-newy),abs(Fz-newz)])/max(abs([Fx,Fy,Fz])));

%     fprintf('Completed force and torque calculation.\n\n');
%     fprintf('The forces are:\nF_x = %f N.\nF_y = %f N.\nF_z = %f N.\n\n',Fx,Fy,Fz);
%     fprintf('The torques are:\ntau_x = %f N.\ntau_y = %f N.\ntau_z = %f N.\n\n\n',taux,tauy,tauz);

%end

F = [Fx,Fy,Fz];
T = [taux,tauy,tauz];

toc;


end



