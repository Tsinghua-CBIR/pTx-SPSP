function aE = ernstAngle(TR,T1)
if nargin == 0
    TR = 500;
    T1 = 1644;
end
if nargin == 1
    T1 = 1939;
end
TR = TR*1e3;
aE = acos(exp(-TR/T1));
aE = aE/pi*180;
end