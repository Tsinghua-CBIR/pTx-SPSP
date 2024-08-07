function lSAR = localSARcom(rf,TR,dt)
load 'data/SarDataUser'
nVOP = size(ZZ,3);

lSAR = zeros(1,nVOP);
for i = 1:nVOP
    Q = ZZ(:,:,i);
    temp = 0;
    for j = 1:size(rf,2)
        temp = temp+rf(:,j)'*Q*rf(:,j);
    end
    lSAR(i) = temp*dt*(1/TR);
end
lSAR = max(abs(lSAR));
%%% When using more than 1000 VOP matrices, we recommend you use a
%%% vectorized manner. As we got only 59 matrices here, we chose to use
%%% loop.
end