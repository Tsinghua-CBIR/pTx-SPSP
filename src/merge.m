function sys = merge(big)
sys = [];
for i = 1:24
    sys(:,i) = sum(big(:,i*12-11:i*12),2);
end
end