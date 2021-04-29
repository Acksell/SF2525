

Ff0_func = @(k) Ff0(k);

Ns=100000;
kspace = linspace(-Ns/100, Ns/100,Ns);
% figure
% plot(arrayfun(@(k) b(abs(k)), kspace)) 
% hold on
% plot(arrayfun(@(k) nu(mu(k)), kspace)) 
% plot(arrayfun(@(k) mu(k), kspace)) 
% figure

Ff_0 = arrayfun(Ff0_func, kspace);

% plot(Ff_0);
% plot(arrayfun(@(k) mu(k), kspace))

f0 = fftshift(ifft(ifftshift(Ff_0)));
plot(imag(f0))

% hold on
% plot(imag(f0),'bo-')
% plot(real(f0),'ro-')

