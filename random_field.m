% Simulate the random field v


x = -100:.01:100;
f = (-10.005:.01:10.005);
Ff0_ = arrayfun(@(f_) Ff0(f_), f);
f0=real(idftransform(Ff0_, x, f))';

% plot(t,real(f0))



% Precompute f_m

m_min=-16;
m_max=4;
M=m_max-m_min+1;


n_min=-10;
n_max=9;
N_bandwidth = n_max-n_min+1;


m_range=m_min:m_max;
n_range=n_min:n_max;

fm=zeros(length(f0),M+1);
for m=1:M
   fm(:,m)=2^(m_range(m)/3)*f0;
end
% figure
% plot(fm)



%% Monte Carlo

denom=0
numer=0



realisations=2000;
for i=1:realisations % realization of vy
%     for x_pos=xlist
    % Precompute gaussians
    gamma=randn(M, N_bandwidth);
    [X,Y,vy]=FractalRandomField(fm, x, gamma, m_range, n_range);

    numer=numer+(vy(38)-vy(1))^4;
    denom=denom+(vy(38)-vy(1))^2;
%     end
end

numer=numer/realisations
denom=(denom/realisations)^2
numer/denom


function [X,Y,vy]=FractalRandomField(fm, x, gamma, m_range, n_range)
    w=0.5; % make sure that w*Tmax < 2^mmax
    dt=0.01;
    t=0:dt:10;
    X=w*t;
    vy=zeros(1,length(X));

    for m=1:length(m_range)
        for n=1:length(n_range)
            arg=2^(-m_range(m))*X - floor(2^(-m_range(m))*X) - n_range(n);
            vy = vy + gamma(m,n)*interp1(x, fm(:,m), arg);
        end
    end

    % Use trapezoidal to integrate v(X) over time.
    Y=zeros(1,length(t));
    for s=2:length(t)
        Y(s)=trapz(0:dt:t(s), vy(1:s));
    end
end
