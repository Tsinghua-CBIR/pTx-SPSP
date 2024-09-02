function cc = SARcomKT(xxx,TR,RFA,dt)
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
cc = RFA*RFA*max(abs(lSAR))*1e12/5; %% Attention, for KT, need to divide nspts, as average the amplitude by time means nspts/(nspts*nspts) = 1/nspts

end