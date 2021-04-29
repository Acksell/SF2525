

Ff0_func = @(k) Ff0(k);


kspace = linspace(1/3, 4/3, 30);
% figure
% plot(arrayfun(@(k) b(abs(k)), kspace)) 
% hold on
% plot(arrayfun(@(k) nu(mu(k)), kspace)) 
% plot(arrayfun(@(k) mu(k), kspace)) 
% figure

Ff_0 = arrayfun(Ff0_func, kspace);

% plot(Ff_0);
% plot(arrayfun(@(k) mu(k), kspace))

f0 = ifft(Ff_0);
plot(abs(f0),'-*')
% plot(real(f0).^2+imag(f0).^2,'*')
% loglog(real(f0), '*')



