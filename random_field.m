% Simulate the random field v


x = -40:.01:40;
f = (-10.005:.01:10.005);
Ff0_ = arrayfun(@(f_) Ff0(f_), f);
f0=real(idftransform(Ff0_, x, f))';

% plot(t,real(f0))

% Precompute f_m

m_min=-40;
m_max=0;
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

X=logspace(log10(2^min(m_range)),log10(2^max(m_range)));
vy_matrix=zeros(length(X),realisations);

for i=1:realisations % realization of vy
    i
%     for x_pos=xlist
    % Precompute gaussians
    gamma_mn=randn(M, N_bandwidth);
    [X,Y,vy]=FractalRandomField(fm, x, X, gamma_mn, m_range, n_range);
    vy_matrix(:,i)=vy;
    numer=numer+(vy(3)-vy(1))^4;
    denom=denom+(vy(3)-vy(1))^2
%     end
end

numer=numer/realisations;
structureFunc=(denom/realisations);
denom=(denom/realisations)^2;
numer/denom

%%
sfunc_list=zeros(length(vy),1);
sfunc_list_real=zeros(length(vy),1);
for pos=2:length(vy)-1
    sfunc=0;
    for i=1:realisations
        sfunc=sfunc+(vy_matrix(pos,i)-vy_matrix(1,i))^2;
    end
    sfunc_list(pos)=sfunc/realisations;
    sfunc_list_real(pos)=2^(8/3)*pi^(5/3)*X(pos)^(2/3)/(sqrt(3)*gamma(5/3));
end
figure
loglog(X,sfunc_list,'x-')
hold on
loglog(X,sfunc_list_real)
legend(["Calculated structure", "Analytic structure"],'Location','nw')
figure
semilogx(X,sfunc_list./sfunc_list_real)
% structureFunc
% realStructureFunc=2^(8/3)*pi^(5/3)*X(3)^(2/3)/(sqrt(3)*gamma(5/3))


function [X,Y,vy]=FractalRandomField(fm, x, X, gamma, m_range, n_range)
    
%     w=0.5; % make sure that w*Tmax < 2^mmax
%     X=w*t;
    vy=zeros(1,length(X));

    for m=1:length(m_range)
        for n=1:length(n_range)
            arg=2^(-m_range(m))*X - floor(2^(-m_range(m))*X) - n_range(n);
            vy = vy + gamma(m,n)*interp1(x, fm(:,m), arg);
        end
    end

    % Use trapezoidal to integrate v(X) over time.
    Y=zeros(1,length(X));
%     for s=2:length(t)
%         Y(s)=trapz(0:dt:t(s), vy(1:s));
%     end
end
