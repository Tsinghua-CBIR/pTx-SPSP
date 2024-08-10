function kpmse = foptKT(myKT)
global A maskMS b1arr posarr nKT gamma b0mapMS phasetrack dt b off trel sysmat myfox wt passband stopband point_P point_S
kp = zeros(3,24);

% for group = 1:8
%     kp(:,group*3-2) = zeros(3,1);
%     kp(:,group*3-1) = -myKT(group*3-2:group*3);
%     kp(:,group*3) = myKT(group*3-2:group*3);
% end

kp = reshape(myKT,[3,24]);


for group = 1:8
    xxx = kp(:,group*3);
    if xxx'*xxx>2500
        kpmse = 99;
        return
    end
end
kpback = [kp zeros(3,1)];
dkval = diff(kpback,1,2);
grad = {};
maxamp = 50e-3;
maxsr = 160;
dt = 10e-6;
nblips = size(dkval,2);
for ind= 1: nblips     
    dkmax = max(abs(dkval(:,ind)));
    %blipmax = design_toptgrad1D(0,0,dkmax,maxamp,maxsr,dt);
    blipmax = design_grad_trapz(dkmax,maxamp,maxsr,dt);
    blips = [blipmax;blipmax;blipmax];
    sf = diag(dkval(:,ind)./ dkmax);
    grad{ind} = sf*blips;  
end

npts = length(grad);
npts0= 3;
for ind=1:npts0:npts
    g0=[];
    for jdx = 1:npts0
        g0= [g0 grad{ind+jdx-1}];
    end
          
%     nt = round((ntimepts - size(g0,2))/npts0); 
    nt = 5;
    if nt>0
        subrfLen(ind:ind+npts0-1) = nt; else
        subrfLen(ind:ind+npts0-1) = 1;
    end
end

phasetrack = subrfLen;
for ind= 1: length(grad)-1
    phasetrack(ind+1) = phasetrack(ind)+ size(grad{ind},2)+ subrfLen(ind+1);
end

[nspa,nchs] = size(b1arr);

TFA = 1;


% A = cell(point_P+point_S,1);
% b = cell(point_P+point_S,1);


pha = posarr*kp;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
off = b0mapMS(maskMS)*kb0(phasetrack);
trel = 1/(gamma/2/pi)*kb0(phasetrack);

for i = 1:point_P
    kernalmat = 1i.* gamma.* dt.* exp(1i* ( pha + off - passband(i)*trel ) );
    for idx = 1:nchs
        iidx0 = (idx-1)*nKT + 1;
        sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
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
        sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
    end
    A((i-1+point_P)*nspa+1:(i+point_P)*nspa,:) = wtS*sysmat;
%     b{i+point_P} = zeros(nspa,1);
end

% A = cell2mat(A);
% b = cell2mat(b);

% 解方程，这里原先直接用的最小二乘，有点草率
 [wt,~,~]= solve_mlstr(A,b,1000,1e-5);
% wt = (A'*A)\(A'*b);

wt1 = wt;
% 1000*(wt1'*wt1)

kpmse = (abs(A*wt1)-abs(b))'*(abs(A*wt1)-abs(b))+4e5*(wt1'*wt1);
% kpmse = (abs(A*wt1)-abs(b))'*(abs(A*wt1)-abs(b))+max(showSAR(wt))/8
% toc;
end