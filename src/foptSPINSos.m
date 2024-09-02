function kpmse = foptSPINSos(spin,wt0)
global A maskMS b1arr posarr nKT gamma b0mapMS b off trel sysmat myfox wt kernalmat passband stopband point_P point_S

Kpp = spin(1);
Ktt = spin(2)/6;
Kphph = spin(3)/30;
alpha = spin(4)/2;
beta = spin(5)/40;

T = 2.88e-3;%sec，total time
dt = 10e-6;%dwell time，sec
len = round(T/dt);
kp = zeros(3,len);
for i = 1:len
    ttt = i*dt;
    Kp = Kpp/(1+exp(alpha*(ttt/T-beta)));
    Ktheta = Ktt*pi*1000*ttt;
    Kphi = Kphph*pi*1000*ttt;
    Kp_temp = [Kp*sin(Ktheta)*cos(Kphi);Kp*sin(Ktheta)*sin(Kphi);Kp*cos(Ktheta)];
    kp(:,i) = Kp_temp;
end

[nspa,nchs] = size(b1arr);

TFA = 1;



% A = cell(point_P+point_S,1);
% b = cell(point_P+point_S,1);


pha = posarr*kp;


for i = 1:point_P
    kernalmat = 1i.* gamma.* dt.* exp(1i* ( pha + off - passband(i)*trel ) );
    for idx = 1:nchs
        iidx0 = (idx-1)*nKT + 1;
        sysmat(:,iidx0:(iidx0+nKT-1)) = (diag(b1arr(:,idx)) * kernalmat);
    end
    A((i-1)*nspa+1:i*nspa,:) = sysmat;
%     tar = (TFA/180)*pi*ones(nspa,1);
%     b{i} = tar;
end

wtS = 8;
for i = 1:point_S
    kernalmat = 1i.* gamma.* dt.* exp(1i* ( pha + off - stopband(i)*trel ) );
    for idx = 1:nchs
        iidx0 = (idx-1)*nKT + 1;
        sysmat(:,iidx0:(iidx0+nKT-1)) = (diag(b1arr(:,idx)) * kernalmat);
    end
    A((i-1+point_P)*nspa+1:(i+point_P)*nspa,:) = wtS*sysmat;
%     b{i+point_P} = zeros(nspa,1);
end

% Solve the reduced problem
 [wt,~,~]= solve_mlstr(A,b,1500,1e-5);

 

% toc;
end