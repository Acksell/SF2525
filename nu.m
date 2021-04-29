function [res]=nu(k)
    res = 2*(max(0, (k-y(0)))^2 + max(0, (k-y(2)))^2 - 2*max(0, k-y(1))^2);
end

function [res]=y(r)
    res=(cos((2-r)*pi/2)+1)/2;
end