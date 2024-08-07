function mse = FAIupdate(xx,vx1)
len = round(length(xx)/2);
x = xx(1:len)+1i*xx(len+1:end);
global A b
step = x-(vx1(1:len)+1i*vx1(len+1:end));
temp = abs(A*x);
mse = (temp-b)'*(temp-b)+1e9*(step'*step);
end