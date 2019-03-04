
% Script to mesh a second magnet and apply a quadric equation to each
% element in order to calculate force
%
% James O'Connell 27th Feb 2019
%
% ToDo: Vectorise this code!

close all;
clear all;
clc;

phi = (1+sqrt(5))/2;
A = 0.01*phi;
verticesA = A*[1,1,1;1,1,-1;1,-1,1;1,-1,-1;-1,1,1;-1,1,-1;-1,-1,1;-1,-1,-1;...
    0,phi,1/phi;0,phi,-1/phi;0,-phi,1/phi;0,-phi,-1/phi;...
    1/phi,0,phi;1/phi,0,-phi;-1/phi,0,phi;-1/phi,0,-phi;...
    phi,1/phi,0;phi,-1/phi,0;-phi,1/phi,0;-phi,-1/phi,0];
thetax = 31.715*pi/180;
R_x = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
verticesA = verticesA*R_x;
magA = [0,0,1.3];
Sa = alphaShape(verticesA,Inf);

verticesB = verticesA + repmat([0,0,max(verticesA(:,3))-min(verticesA(:,3))+0.01],length(verticesA),1);
magB = [0,0,-1.3];
Sb = alphaShape(verticesB,Inf);

% Calculate the actual forces using current force method
[Factual,Tactual] = polyhedronForce(verticesA,verticesB,magA,magB,mean(verticesB),1e-30,0.45);
Factual;

% Now let's estimate the force on magnet B
Fac = minConvexHull(verticesB);
[Ver,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
[Fac,~] = triangulateFaces(Fac);
[Ver,Fac] = subdivideMesh(Ver,Fac,4);
norms = meshFaceNormals(Ver,Fac);
Areas = meshFaceAreas(Ver,Fac);

Frough = [0,0,0];
Fest = [0,0,0];
for k = 1:size(Fac,1)
    MdotN = myDot(norms(k,:),magB);
    n = norms(k,:)/norm(norms(k,:));
    thetay = atan2(n(1),-n(3));
    thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
    Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
    Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
    Rn = Rx*Ry;
    R = Rn^-1;
    
    facepts = Ver(Fac(k,:),:);
    faceptsp = [facepts;mean(facepts(1:2,:));mean(facepts(2:3,:));mean(facepts([1,3],:))];
    fieldsolns = polyhedronField(verticesA,magA,faceptsp);
    fieldz = fieldsolns(:,3);
    facerotp = faceptsp*R;
    facerot = facepts*R;
    fieldrot = fieldsolns*R;
    xi = facerotp(:,1);
    yi = facerotp(:,2);
    Amat = [xi.^2,xi,yi.^2,yi,xi.*yi,repmat([1],6,1)];
    Xmat = Amat\fieldrot;
    
    mypts = mean(facerotp);
    fieldestimate = [mypts(1)^2,mypts(1),mypts(2)^2,mypts(2),mypts(1)*mypts(2),1]*Xmat;
%     MdotN*(fieldestimate*Areas(k))*Rn / (4*pi*10^-7)
    Frough = Frough + MdotN*(fieldestimate*Areas(k))*Rn / (4*pi*10^-7);
    
    FF = 0;
%     facerot = sortrows(facerot);
    % First sub-triangle:
    m = [(facerot(2,2)-facerot(1,2))/(facerot(2,1)-facerot(1,1));(facerot(3,2)-facerot(1,2))/(facerot(3,1)-facerot(1,1))];
%     m = sort(m);
    c = facerot(1,2)-m*facerot(1,1);
    x = facerot(1:2,1);
    for i = 1:2
        for j = 1:2
%             FF = FF + (-1)^(i+j)*(1/24*x(j)*(4*c(i)*(2*Xmat(3,:)*c(i)^2+3*c(i)*Xmat(4,:)+6*Xmat(6,:))+6*(2*Xmat(2,:)*c(i)+2*c(i)*Xmat(4,:)*m(i)+2*Xmat(6,:)*m(i)+c(i)^2*(Xmat(5,:)*2*Xmat(3,:)*m(i)))*x(j)+4*(2*Xmat(1,:)*c(i)+m(i)*(2*Xmat(2,:)+Xmat(4,:)*m(i)+2*c(i)*(Xmat(5,:)+Xmat(3,:)*m(i))))*x(j)^2+m(i)*(6*Xmat(1,:)+m(i)*(3*Xmat(5,:)+2*Xmat(3,:)*m(i)))*x(j)^3));
            FF = FF + (-1)^(i+j)*1/24*x(j)*(4*c(i)*(2*Xmat(3,:)*c(i)^2+3*c(i)*Xmat(4,:)+6*Xmat(6,:))+6*(2*Xmat(2,:)*c(i)+2*c(i)*Xmat(4,:)*m(i)+2*Xmat(6,:)*m(i)+c(i)^2*(Xmat(5,:)+2*Xmat(3,:)*m(i)))*x(j)+4*(2*Xmat(1,:)*c(i)+m(i)*(2*Xmat(2,:)+Xmat(4,:)*m(i)+2*c(i)*(Xmat(5,:)+Xmat(3,:)*m(i))))*x(j)^2+m(i)*(6*Xmat(1,:)+m(i)*(3*Xmat(5,:)+2*Xmat(3,:)*m(i)))*x(j)^3);
%             FF = FF + (-1)^(i+j)*(x(j)^4/4*(Xmat(1,:)*m(i)+Xmat(3,:)*m(i)^3/3)+x(j)^3/3*(Xmat(1,:)*c(i)+Xmat(2,:)*m(i)+Xmat(3,:)*m(i)^2*c(i)+Xmat(4,:)*m(i)^2/2+Xmat(5,:)*m(i)^2/2)+x(j)^2/2*(Xmat(2,:)*c(i)+Xmat(3,:)*m(i)*c(i)+Xmat(4,:)*m(i)*c(i)+Xmat(5,:)*m(i)*c(i)+Xmat(6,:)*m(i))+x(j)*(Xmat(3,:)*c(i)^3/3+Xmat(4,:)*c(i)^2/2+Xmat(5,:)*c(i)^2/2+Xmat(6,:)*c(i)));
        end
    end
    
%     MdotN*FF*Rn / (4*pi*10^-7)
    % Second sub-triangle:
    m = [(facerot(3,2)-facerot(2,2))/(facerot(3,1)-facerot(2,1));(facerot(3,2)-facerot(1,2))/(facerot(3,1)-facerot(1,1))];
%     m = sort(m);
    c = facerot(3,2)-m*facerot(3,1);
    x = facerot(2:3,1);
    for i = 1:2
        for j = 1:2
%             FF = FF + (-1)^(i+j)*(1/24*x(j)*(4*c(i)*(2*Xmat(3,:)*c(i)^2+3*c(i)*Xmat(4,:)+6*Xmat(6,:))+6*(2*Xmat(2,:)*c(i)+2*c(i)*Xmat(4,:)*m(i)+2*Xmat(6,:)*m(i)+c(i)^2*(Xmat(5,:)*2*Xmat(3,:)*m(i)))*x(j)+4*(2*Xmat(1,:)*c(i)+m(i)*(2*Xmat(2,:)+Xmat(4,:)*m(i)+2*c(i)*(Xmat(5,:)+Xmat(3,:)*m(i))))*x(j)^2+m(i)*(6*Xmat(1,:)+m(i)*(3*Xmat(5,:)+2*Xmat(3,:)*m(i)))*x(j)^3));
            FF = FF + (-1)^(i+j)*1/24*x(j)*(4*c(i)*(2*Xmat(3,:)*c(i)^2+3*c(i)*Xmat(4,:)+6*Xmat(6,:))+6*(2*Xmat(2,:)*c(i)+2*c(i)*Xmat(4,:)*m(i)+2*Xmat(6,:)*m(i)+c(i)^2*(Xmat(5,:)+2*Xmat(3,:)*m(i)))*x(j)+4*(2*Xmat(1,:)*c(i)+m(i)*(2*Xmat(2,:)+Xmat(4,:)*m(i)+2*c(i)*(Xmat(5,:)+Xmat(3,:)*m(i))))*x(j)^2+m(i)*(6*Xmat(1,:)+m(i)*(3*Xmat(5,:)+2*Xmat(3,:)*m(i)))*x(j)^3);
%             FF = FF + (-1)^(i+j)*(x(j)^4/4*(Xmat(1,:)*m(i)+Xmat(3,:)*m(i)^3/3)+x(j)^3/3*(Xmat(1,:)*c(i)+Xmat(2,:)*m(i)+Xmat(3,:)*m(i)^2*c(i)+Xmat(4,:)*m(i)^2/2+Xmat(5,:)*m(i)^2/2)+x(j)^2/2*(Xmat(2,:)*c(i)+Xmat(3,:)*m(i)*c(i)+Xmat(4,:)*m(i)*c(i)+Xmat(5,:)*m(i)*c(i)+Xmat(6,:)*m(i))+x(j)*(Xmat(3,:)*c(i)^3/3+Xmat(4,:)*c(i)^2/2+Xmat(5,:)*c(i)^2/2+Xmat(6,:)*c(i)));
        end
    end
%     MdotN*FF*Rn / (4*pi*10^-7)
    Fest = Fest + MdotN*FF*Rn / (4*pi*10^-7);
end
Frough = Frough;
Fest = Fest;
toc

figure;
drawMesh(Ver,Fac);

Fact = [-0.0001,-0.0016,266.2348]; % 1cm gap
% Fact = [0.0015   -0.0103  547.7106]; % 1mm gap
(Factual(end,3)-Fact(3))/Fact(3)*100
(Fest(3)-Fact(3))/Fact(3)*100



