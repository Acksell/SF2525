% Simulate the random field v


x = -20:.01:20;
f = (-10.005:.01:10.005);
Ff0_ = arrayfun(@(f_) Ff0(f_), f);
f0=real(idftransform(Ff0_, x, f))';

% plot(t,real(f0))



% Precompute f_m

m_min=-8;
m_max=3;
M=m_max-m_min+1;


n_min=-10;
n_max=9;
N_bandwidth = n_max-n_min+1;

fm=zeros(length(f0),M+1);
for m=1:M+1
   fm(:,m)=2^(m/3)*f0;
end

% plot(fm)


% Precompute gaussians
gamma=randn(M, N_bandwidth);
m_range=m_min:m_max;
n_range=n_min:n_max;

minarg=zeros(N_bandwidth,M);
maxarg=zeros(N_bandwidth,M);

w=10;
dt=0.01;
t=0:dt:10;
X=w.*t;
vy=zeros(1,length(X));

for m=1:length(m_range)
    for n=1:length(n_range)
        arg=2^(-m_range(m))*X-floor(2^(-m_range(m))*X) - n_range(n);
        vy = vy + gamma(m,n)*interp1(x, fm(:,m), arg);
    end
end

% Use trapezoidal to integrate v(X) over time.
Y=zeros(1,length(t));
for s=2:length(t)
    Y(s)=trapz(0:dt:t(s), vy(1:s));
end

% plotting

plot(vy)
% plot(Y)
% plot(X,Y);
%%
for m=1:length(m_range)
    arginterp=2^(-m_range(m))*X-floor(2^(-m_range(m))*X)-n_range(n);
    figure 
    plot(arginterp)
%     figure
%     plot(arginterp, interp1(x, fm(:,m), arginterp),'o')
%     hold on
%     plot(x, fm(:,m))
end





