%% Monte Carlo simulation of a Kou jump-diffusion process

% Parameters
clear all 

% Market Parameter
T = 1; % maturity
S0 = 1; % spot price
K = 1.1; % strike price
q = 0.02; % dividend rate
r = 0.05; % risk-free interest rate
muS = 0.2; sigma = 0.3; % model parameters of the diffusion part (GBM)

% Model Parameter 
npaths = 200000; % number of paths
n = 10^6;
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
eta1 = 4; % parameter of the exponential distribution controlling the upward jumps
eta2 = 3; % parameter of the exponential distribution controlling the downward jumps
p = 0.4; % probability of an upward jump
lambda = 0.5; 
binw = 0.1; % bin width
xmax = 2; % truncation
deltax = 0.01; % grid step

%% Monte Carlo

% Compute the increments of the ABM
dW = (muS-0.5*sigma^2)*dt + sigma*sqrt(dt)*randn(nsteps,npaths);

%%Compute the PDF and sample the distribution
x = -xmax:deltax:xmax; % grid
fX = p*eta1*exp(-eta1*x).*(x>=0) + (1-p)*eta2*exp(eta2*x).*(x<0); % PDF

% Compute the increments of the double exponential distribution
U = rand(nsteps,npaths); % standard uniform random variable
X = -1/eta1*log((1-U)/p).*(U>=1-p)+1/eta2*log(U/(1-p)).*(U<1-p); % bilateral exp. r.v.

% Compute the increments of the NCPP
dN = poissrnd(lambda*dt,[nsteps,npaths]); % Find our Poisson random number
dJ = dN.*X;

% Sum the increments of the ABM and the NCPP
dX = dW + dJ;

% Accumulate the increments
X = [zeros(1,npaths); cumsum(dX)];

S = S0*exp(X);

VcMCb = exp(-r*T)*mean(max(S(end,:)-K,0),2)
VpMCb = exp(-r*T)*mean(max(K-S(end,:),0),2)

%% Expected mean and sample path

figure(1)
close all
EX = (muS-0.5*sigma^2+lambda*(p/eta1-(1-p)/eta2))*t; % expected path
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('S')
ylim([-1,1.2]);
title('Paths of a kou jump-diffusion process X = \mut + \sigmaW(t) + \Sigma_{i=1}^{N(t)} Z_i')
print('-dpng','mjdpaths.png')

%% Variance = mean square deviation (Analytic expression for variance)

figure(2)
plot(t,(sigma^2+2*lambda*(p/eta1^2+(1-p)/eta2^2))*t,'b', t, var(X,0,2),'r',t,mean((X-transpose(EX)).^2,2),'g')
legend('Theory', "Sampled1", "Sampled2")
xlabel('t')
ylabel('Var(X) = E((X-E(X))^2)')
title('Kou jump-diffusion process X = \mut + \sigmaW(t) + \Sigma_{i=1}^{N(t)} Z_i: variance')

%% Probability density function at different times

figure(3)

[h,x] = hist(X(40,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,1)
bar(x,f)
ylabel('f_X(x,0.2)')
xlim([-1,1.2]);
ylim([0,3]);
title('Probability density function of a Kou jump-diffusion process at different times')

[h,x] = hist(X(100,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,2)
bar(x,f);
xlim([-1,1.2]);
ylim([0,3]);
ylabel('f_X(x,0.5)')

[h,x] = hist(X(end,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,3);
bar(x,f);
xlim([-1,1.2]);
ylim([0,3]);
xlabel('x')
ylabel('f_X(x,1)')

%% Monte Carlo - Blocks

tic
nblocks = 20; % number of blocks
newsample = 10000; % number of samples per block

VcMCb = zeros(nblocks,1);
VpMCb = zeros(nblocks,1);
for i = 1:nblocks
    % Compute the increments of the ABM
    dW = (muS-0.5*sigma^2)*dt + sigma*sqrt(dt)*randn(nsteps,newsample);
    
    % Compute the increments of the double exponential distribution
    U = rand(nsteps,newsample); % standard uniform random variable
    dZ = -1/eta1*log((1-U)/p).*(U>=1-p)+1/eta2*log(U/(1-p)).*(U<1-p); % bilateral exp. r.v.
    
    % Compute the increments of the NCPP
    dN = poissrnd(lambda*dt,[nsteps,newsample]); % Find our Poisson random number
    
    % Sum the increments of the ABM and the NCPP
    dX = dW + dZ.*dN;
    
    % Accumulate the increments
    X = [zeros(1,newsample); cumsum(dX)];

    % Stock price
    S = S0*exp(X);
    
    VcMCb(i) = exp(-r*T)*mean(max(S(end,:)-K,0),2);
    VpMCb(i) = exp(-r*T)*mean(max(K-S(end,:),0),2);
end     

VcMC = mean(VcMCb)
VpMC = mean(VpMCb)
cputime_MC = toc;
scMC = sqrt(var(VcMCb)/nblocks);
spMC = sqrt(var(VpMCb)/nblocks);

%% Analyticl solution

muRN = r-q-0.5*sigma^2; % drift coefficient of the arithmetic Brownian motion
d2 = (log(S0/K)+muRN*T)/(sigma*sqrt(T));
d1 = d2 + sigma*sqrt(T);
Vca = S0*exp(-q*T)*cdf('Normal',d1,0,1) - K*exp(-r*T)*cdf('Normal',d2,0,1)
Vpa = K*exp(-r*T)*cdf('Normal',-d2,0,1) - S0*exp(-q*T)*cdf('Normal',-d1,0,1)
