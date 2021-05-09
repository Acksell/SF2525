function [res]=Fphi(k)
    res=1i*sign(-k)*exp(1i*pi*k)*b(abs(k));
end