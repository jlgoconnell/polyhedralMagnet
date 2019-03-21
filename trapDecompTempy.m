
% Function to decompose a polygon in the XY plane into a set of trapezia
% with the parallel sides parallel to the Y-axis.
%
% James O'Connell 20th March 2019


function trapVertices = trapDecomp(vertices)

close all;

[~,I] = sort(vertices(:,1));
vertices = vertices(I,:);

conlist = convhull(vertices);
vertices = vertices(conlist(1:end-1),:);

s = alphaShape(vertices,inf);
subplot(1,2,1);
cc = convhull(s.Points);
fill(s.Points(cc,1),s.Points(cc,2),'m');
xlim([min(vertices(:,1))-0.25,max(vertices(:,1))+0.25]);
ylim([min(vertices(:,2))-0.25,max(vertices(:,2))+0.25]);
axis square;
hold on;
% for j = 1:length(conlist)-1
%     text(vertices(j,1),vertices(j,2),num2str(j));
% end

[x,I] = sort(vertices(:,1));

I = [I;I(end-1:-1:1)];
ms = [];%NaN(length(I)-1,1);
% cs = ms;
js = ms;

for j = 1:max(I)
    Iind = length(ms)+1;%find(ms==NaN,1);
    js(find(I(Iind:end)==j,1)+Iind-1:find(I(Iind:end)==mod(j,max(I))+1,1)+Iind-2,1) = j;
    ms(find(I(Iind:end)==j,1)+Iind-1:find(I(Iind:end)==mod(j,max(I))+1,1)+Iind-2,1) = (vertices(mod(j,max(I))+1,2)-vertices(j,2))/(vertices(mod(j,max(I))+1,1)-vertices(j,1));
end
cs = vertices(js,2)-ms.*vertices(js,1);
ms = [ms(1:length(ms)/2),ms(end:-1:length(ms)/2+1)]';
cs = [cs(1:length(cs)/2),cs(end:-1:length(cs)/2+1)]';
% js = [js(1:length(js)/2),js(end:-1:length(js)/2+1)]';
% x = x(js);

ms = [ms,ms(:,end)];
cs = [cs,cs(:,end)];
xs = repmat(x',2,1);

ind = ~sum(isinf(ms),1);
ms = ms(:,ind);
cs = cs(:,ind);
xs = xs(:,ind);


trapVertices = ms.*xs+cs;
ind = ~prod(diff([trapVertices';NaN,NaN])==0,2)
trapVertices = trapVertices(:,ind);
xs = xs(:,ind);

subplot(1,2,2);
hold on;
these = trapVertices(:);
xx = xs(:);
for j = 1:2:length(these)-2
    fill(xx(j:j+3),these([j,j+1,j+3,j+2]),rand(1,3))
end
xlim([min(vertices(:,1))-0.25,max(vertices(:,1))+0.25]);
ylim([min(vertices(:,2))-0.25,max(vertices(:,2))+0.25]);
axis square;

end







