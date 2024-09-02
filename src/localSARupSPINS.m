function x1 = localSARupSPINS(x0,TR,RFA,dt)
lS = SARcomSPINS(x0,TR,RFA,dt);
if lS>9 && lS<10.000000001
%%%If it happens to fall within this range, there is almost no room for optimization and you can exit directly
x1 = x0;
return
end
len = round(length(x0)/2);
xx = x0(1:len)+1i*x0(len+1:end);
global A b
if lS<10
   [x1,~,~]= solve_mlstr(A,b,1000/10*(lS+1),1e-5,xx);
   x1 = [real(x1);imag(x1)];
else
    penalty = lS/10;
    x2 = x0;
    while (SARcomSPINS(x2,TR,RFA,dt)>10)
        penalty = penalty+lS/10;
        xx = x2(1:len)+1i*x2(len+1:end);
        [xx,~,~]= solve_mlstr(A,b,1000/10*(lS+penalty),1e-5,xx);
        x2 = [real(xx);imag(xx)];
    end
    
    xx = x0(1:len)+1i*x0(len+1:end);
    
    load 'data/SarDataUser' ZZ
    nVOP = size(ZZ,3);
    nchs = size(ZZ,1);
    rf = reshape(xx,[],nchs);
    rf = rf.';
    xxx = sum(rf,2);
    xxx0 = xxx;
    lSAR = zeros(1,nVOP);
    for i = 1:nVOP
        Q = ZZ(:,:,i);
        
        lSAR(i) = xxx'*Q*xxx;
    end
    [~,index] = max(lSAR);
    Q = dt*(1/TR)*RFA*RFA*1e12*ZZ(:,:,index);
    miu = 0;
    R = eye(nchs);
    F = 10;
    xxx0 = xxx;
    for i = 1:10
        xn = inv(R+miu*Q)*(xxx0);
        mse = xn'*Q*xn;
        xnew = inv(R+(miu+1e-10)*Q)*xxx0;
        delta = mse-xnew'*Q*xnew
        if delta<1e-10 || mse<F
            break
        end
        xxx = xn;
        miu = miu+(mse-F)/delta;
    end
    rf = rf.*repelem(xxx./xxx0,1,24);
    x3 = reshape(rf',[],1);
    x3 = [real(x3);imag(x3)];
    xx2 = x2(1:len)+1i*x2(len+1:end);
    xx3 = x3(1:len)+1i*x3(len+1:end);
    if (A*xx2-b)'*(A*xx2-b)<(A*xx3-b)'*(A*xx3-b)
        x1 = [real(x2);imag(x2)];
    else
        x1 = [real(x3);imag(x3)];
    end
end
end