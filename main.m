
% Ns=101;
% Note: Ns/100 is simply found by trial and error.
% Note also: Ff0 has support on the interval [1/3, 4/3]

N = 1024;
fs = 10;
df = fs / N;
f = linspace(-N*df/2, N*df/2, N);
wspace = f;

% wspace = linspace(-10, 10, N);
% wspace=kspace*2*pi;

% figure
% plot(arrayfun(@(k) b(abs(k)), kspace)) 
% hold on
% plot(arrayfun(@(k) nu(mu(k)), kspace)) 
% plot(wspace, arrayfun(@(w) mu(w/(2*pi)), wspace)) 
% figure

% Generate the fourier-transformed values of f0.
Ff_0 = arrayfun(@(w) Ff0(w), wspace);

Ff_1 = arrayfun(@Ff1, wspace);
f1=fftshift(ifft(ifftshift(Ff_1)));


% fftshift is needed to reorder the values before use of ifft.
% ifftshift shifts them back.


% NOTE: maybe some weirdness becayse ifft expects frequency and not k?

f0 = fftshift(ifft(ifftshift(Ff_0)));
% Problem: Imag looks like the correct f0 function, but in reality we want
% f0 to be a real function and for the imaginary part to simply be rounding errors.
figure
plot(real(f0))
title("Real part")

figure
plot(imag(f0))
title("Imaginary part")


% figure
% plot(arrayfun(@nu, linspace(-1,2,100)))

figure
F_phi = arrayfun(@Fphi, wspace);
plot(wspace, abs(F_phi))
hold on

title("abs(F\_phi)") % Spectrum of meyer wavelet, matches the wikipedia article.


function [res]=Ff1(k)
    res=exp(-k*k);
end

% function [res]=Fphi(k)
%     res=sign(k)*exp(1i*pi*k)*b(abs(k))/sqrt(2*pi);
% end
% 
% function [res]=Ff0(k)
%     res=sqrt(E(k))*Fphi(k);
% end
% 
% function [res]=nu(k)
%     res = 2*(max(0, (k-y(0)))^2 + max(0, (k-y(2)))^2 - 2*max(0, k-y(1))^2);
% end
% 
% function [res]=y(r)
%     res=(cos((2-r)*pi/2)+1)/2;
% end
% 
% function [res]=b(k)
%     res=sin(nu(mu(k))*pi/2);
% end
% 
% function [res]=E(k)
%     res=abs(k)^(-5/3);
% end
% 
% function [res]=mu(k)
%     absk=abs(k);
%     if absk<1/3
%         res=0;
%     elseif absk<2/3 && 1/3<=absk
%         res=3*(absk-1/3);
%     elseif absk<4/3 && 2/3 <= absk
%         res=1 + 3/2*(2/3-absk);
%     else
%         res=0;
%     end
% end
