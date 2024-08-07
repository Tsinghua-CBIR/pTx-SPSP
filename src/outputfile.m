function [] =  outputfile(rf,grad);
[Nc,Nt] = size(rf);
rf = reshape(rf',[Nc*Nt,1]);
temp = grad;
temp(1,:) = grad(1,:);
temp(2,:) = grad(2,:);
temp(3,:) = grad(3,:);
grad = temp';

save_pTXRFPulse_toINI(grad,rf,[]);