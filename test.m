% inverse fourier
t=(-2):.01:2;
f=(-5):.1:5;
% f=2*pi*f

g=exp(-f.*f);

X=idftransform(g, t, f);

figure;
% stem(t,real(X)) % numerical approximation
plot(t, real(X),'o')
hold on
plot(t, exp((-t.^2/4))/sqrt(2)) % analytical solution

%%
% inverse fourier
t = -5:.01:5;
f = (-10.005:.01:10.005);
Ff0_ = arrayfun(@(f_) Ff0(f_), f);
Fphi_ = arrayfun(@(f_) Fphi(f_), f);
X=idftransform(Ff0_, t, f);

figure;
plot(t, real(X))
title("real")
figure
plot(t, imag(X))
title("imag")


% figure
% plot(imag(Ff0_))
% hold on
% plot(real(Ff0_))
% legend(["imaginary","real"])
% title("Ff0")


% figure
% plot(f,imag(Fphi_))
% hold on
% plot(f,real(Fphi_))
% legend(["imaginary","real"])
% title("Fphi")
