function kpmse = funKT(enlargeX)
global A maskMS b1arr posarr nKT gamma b0mapMS phasetrack dt b off trel sysmat myfox wt passband stopband point_P point_S

len = round(size(enlargeX,1)/2);
wt1 = enlargeX(1:len)+1i*enlargeX(len+1:2*len);
kpmse = 50*(abs(A*wt1)-abs(b))'*(abs(A*wt1)-abs(b));
% kpmse = (abs(A*wt1)-abs(b))'*(abs(A*wt1)-abs(b))+max(showSAR(wt))/8
% toc;
end