


function Force = polyhedronForceQuadric(verticesA,verticesB,magA,magB,nSubdivide)

Force = [0,0,0];

% Get face and vertex information of the polyhedra
FacA = minConvexHull(verticesA);
FacB = minConvexHull(verticesB);
[VerB,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
norms = meshFaceNormals(VerB,FacB);
MdotN = 1/(pi*4e-7)*dot(repmat(magB,size(norms,1),1)',norms')';
FacB = FacB(abs(MdotN)>eps);
norms = meshFaceNormals(VerB,FacB);

% Take each facet and rotate to be parallel to XY axis
for i = 1:length(FacB)
    
    Forcepoly = zeros(1,3);
    
    faceverts = VerB(FacB{i},:);
    n = norms(i,:)/norm(norms(i,:));
    thetay = atan2(n(1),-n(3));
    thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
    Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
    Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
    Rn = Rx*Ry;
    R = Rn^-1;
    faceverts = faceverts*R;
    z = faceverts(1,3);
    vertTrap = trapDecomp(faceverts(:,1:2));
    
    % Take each trapezium and calculate the force on it:
    for j = 1:2:length(vertTrap)-3
        
        myverts = vertTrap(j:j+3,:);
        [~,ia,~] = unique(myverts,'rows');
        if length(ia) == 3
            f = ia';
        else
            f = [1,2,3;2,3,4];
        end
        % Subdivide mesh:
        if nSubdivide >= 2
            [v2,f2] = subdivideMesh(myverts,f,nSubdivide);
        else
            v2 = myverts;
            f2 = f;
        end
        % Calculate the 6 points needed for each triangle:
        [v3,f3] = subdivideMesh(v2,f2,2);
        f4 = zeros(length(f3)/4,6);
        for k = 1:length(f3)/4
            f4(k,1:6) = unique(f3(4*k-3:4*k,:))';
        end
        fieldpts = [v3(f4,:),z*ones(numel(f4),1)];
        % Calculate the field at each point:
        verticesAtemp = verticesA*R;
        Bfield = polyhedronField(verticesAtemp,magA,fieldpts,FacA);
        
        % Set up sparse matrix:
        A = sparse(6*size(f2,1),6*size(f2,1));
        for k = 1:size(f4,1)
            x = v3(f4(k,:),1);
            y = v3(f4(k,:),2);
            A(6*k-5:6*k,6*k-5:6*k) = [x.^2,x,y.^2,y,x.*y,ones(6,1)];
        end
        
        % Solve for coefficients:
        coeffs = A\Bfield;
        A = repelem(coeffs(1:6:end,:),2,2);
        B = repelem(coeffs(2:6:end,:),2,2);
        C = repelem(coeffs(3:6:end,:),2,2);
        D = repelem(coeffs(4:6:end,:),2,2);
        E = repelem(coeffs(5:6:end,:),2,2);
        F = repelem(coeffs(6:6:end,:),2,2);
        coeffsrep = repelem(coeffs,2,2);
        
        % Set up m and c matrices:
        m = [(v2(f2(:,1),2)-v2(f2(:,2),2))./(v2(f2(:,1),1)-v2(f2(:,2),1)),...
            (v2(f2(:,1),2)-v2(f2(:,3),2))./(v2(f2(:,1),1)-v2(f2(:,3),1)),...
            (v2(f2(:,3),2)-v2(f2(:,2),2))./(v2(f2(:,3),1)-v2(f2(:,2),1))];
        mm = m';
        [~,Imax] = max(abs(mm));
        Imax = 4-Imax'+(0:length(Imax)-1)'*size(f2,2);
        ff = f2';
        indices = ff(Imax);
        mm = reshape(mm(abs(mm)~=max(abs(mm))),2,numel(mm)/3);
        m = sort(mm);
        m = m(:);
        c = repelem(v2(indices,2),2,1)-m.*repelem(v2(indices,1),2,1);
        m = repmat(m,1,6);
        c = repmat(c,1,6);
        xs = reshape(v2(f2',1),size(f2'));
        x = [min(xs)',max(xs)'];
        x = repelem(x,2,1);
        x = repmat(x,1,3);
        I = repmat([1,-1;-1,1],size(f2,1),3);
        
        % Estimate force from trapezial plate:
        ftrap = I.*(x.^4/4.*(A.*m+c.*m.^3/3)+x.^3/3.*(A.*c+B.*m+C.*m.^2.*c+D.*m.^2/2+E.*m.^2/2)+x.^2/2.*(B.*c+C.*m.*c+D.*m.*c+E.*m.*c+F.*m)+x.*(C.*c.^3/3+D.*c.^2/2+E.*c.^2/2+F.*c));
        ftrapsum = sum(reshape(ftrap,2*size(ftrap,1),3))
        
        Forcepoly = Forcepoly + ftrapsum;
        
    end
    Forcepoly = Forcepoly*Rn;
    
end







end