function cc = SARcomSPINS(xxx,TR,RFA,dt)
len = round(length(xxx)/2);
xx = xxx(1:len)+1i*xxx(len+1:end);
load 'data/SarDataUser'
nVOP = size(ZZ,3);
nchs = size(ZZ,1);
rf = reshape(xx,[],nchs);
rf = rf.';

lSAR = zeros(1,nVOP);
for i = 1:nVOP
    Q = ZZ(:,:,i);
    temp = 0;
    for j = 1:size(rf,2)
        temp = temp+rf(:,j)'*Q*rf(:,j);
    end
    lSAR(i) = temp*dt*(1/TR);
end
cc = RFA*RFA*max(abs(lSAR))*1e12;

%%% When using more than 1000 VOP matrices, we recommend you use a
%%% vectorized manner. As we got only 59 matrices here, we chose to use
%%% loop.
end