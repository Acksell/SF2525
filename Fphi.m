function [res]=Fphi(k)
    res=sign(k)*exp(1i*pi*k)*b(abs(k))/sqrt(2*pi);
end