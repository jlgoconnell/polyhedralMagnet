
% Script to create some cube magnets in attempt to find a linear magnetic
% spring

% James O'Connell 12th June 2019

function r2 = linearspring_20190612(X)

offx = X(1);
offz = X(2);

magsize = 0.01;

offx = offx + magsize;
% offx = 0.02+magsize; % The second term is to make sure the top and offset magnets don't intersect
% offz = -0.04;
dnominal = 0.2;

magtop = magnetdefine('type','cuboid','dim',[magsize,magsize,magsize],'magn',1,'magdir',[0,0,1]);
magbot = magnetdefine('type','cuboid','dim',[magsize,magsize,magsize],'magn',1,'magdir',[0,0,1]);
magoft = magnetdefine('type','cuboid','dim',[magsize,magsize,magsize],'magn',1,'magdir',[1,0,0]);
magofb = magnetdefine('type','cuboid','dim',[magsize,magsize,magsize],'magn',1,'magdir',[1,0,0]);
magflo = magnetdefine('type','cuboid','dim',[magsize,magsize,magsize],'magn',1,'magdir',[0,0,-1]);

N = 500;
offtop = repmat([0;0;dnominal],1,N);
offbot = repmat([0;0;-dnominal],1,N);
offoft = repmat([offx;0;dnominal+offz],1,N);
offofb = repmat([offx;0;-dnominal-offz],1,N);

n = 4;
displ = linspace(-dnominal+n*magsize,dnominal-n*magsize,N+2);
displ = displ(2:end-1);
displtop = offtop + [0;0;1]*displ;
displbot = offbot + [0;0;1]*displ;
disploft = offoft + [0;0;1]*displ;
displofb = offofb + [0;0;1]*displ;

Ftop = magnetforces(magtop,magflo,displtop);
Fbot = magnetforces(magbot,magflo,displbot);
Foft = magnetforces(magoft,magflo,disploft)*2;
Fofb = magnetforces(magofb,magflo,displofb)*2;

Ftotal = Ftop + Fbot + Foft + Fofb;
Fz = Ftotal(3,:);

% figure;
plot(displ,Fz);
grid on;
hold on;
% 
p = polyfit(displ,Fz,1);
est = p(1)*displ+p(2);
plot(displ,est);
hold off;

r = corrcoef(Fz,est);
r2 = 1-r(2,1)^2
% r2 = max(abs(diff(Fz)/(displ(2)-displ(1))))

end


