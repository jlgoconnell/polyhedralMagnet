
% Script to mesh with trapDecomp.m
%
% James O'Connell 2nd April 2019

close all;
clear;
clc;

vertices = [0.0182,0.0083,0.0233; ...
    0.0120,-0.0108,0.0233; ...
    -0.0080,-0.0108,0.0233; ...
    -0.0142,0.0083,0.0233; ...
    0.0020,0.0200,0.0233];

vertTrap = trapDecomp(vertices(:,1:2));
z = vertices(1,3);

nSubdiv = 3;

for j = 1:2:length(vertTrap)-3 % For each trapezium
    
    myverts = vertTrap(j:j+3,:);
    
    [a,ia,~] = unique(myverts,'rows');
    if length(a) == 3
        f = ia';
    else
        f = [1,2,3;2,3,4];
    end
    
    [v2,f2] = subdivideMesh(myverts,f,nSubdiv); % This is the mesh
    
    [v3,f3] = subdivideMesh(v2,f2,2);
    f4 = []; % This is the list of the 6 points for each triangle
    for k = 1:length(f3)/4
        f4(k,1:6) = unique(f3(4*k-3:4*k,:))';
    end
    
    fieldpts = [v3(f4,:),z*ones(numel(f4),1)];
    
end