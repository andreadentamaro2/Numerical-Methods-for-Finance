%% Monte Carlo simulation of a Kou jump-diffusion process

% Parameters
clear all % clear all variables from memory
npaths = 20000; % number of paths
n = 10^6;
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
eta1 = 4; % parameter of the exponential distribution controlling the upward jumps
eta2 = 3; % parameter of the exponential distribution controlling the downward jumps
p = 0.4; % probability of an upward jump
lambda = 0.5
S0 = 1; 
muS = 0.2; sigma = 0.3; % model parameters of the diffusion part (GBM)
xmax = 2; % truncation
deltax = 0.01; % grid step
binw = 0.1; % bin width

%% Compute the PDF and sample the distribution 
x = -xmax:deltax:xmax; % grid
fX = p*eta1*exp(-eta1*x).*(x>=0) + (1-p)*eta2*exp(eta2*x).*(x<0); % PDF only simulate asymetric bilateral probability density function

U = rand(1,n); % standard uniform random variable
Z = -1/eta1*log((1-U)/p).*(U>=1-p)+1/eta2*log(U/(1-p)).*(U<1-p); % bilateral exp. r.v.


%% Monte Carlo

% Compute the increments of the ABM
dW = (muS-0.5*sigma^2)*dt + sigma*sqrt(dt)*randn(nsteps,npaths); 

% Compute the increments of the double-sided exponential distribution
U = rand(nsteps,npaths); % standard uniform random variable
dZ = -1/eta1*log((1-U)/p).*(U>=1-p)+1/eta2*log(U/(1-p)).*(U<1-p);% bilateral exp. r.v.

% Compute the increments of the NCPP
dN = poissrnd(lambda*dt,[nsteps,npaths]); % Find our Poisson random number

% Sum the increments of the ABM and the NCPP
dX = dW + dZ.*dN;

% Accumulate the increments
X = [zeros(1,npaths); cumsum(dX)];

%% Plot
close all;
figure(1)
x2 = -xmax:binw:xmax; % bin edges
histogram(Z,x2,'normalization','pdf');
hold on
plot(x,fX,'linewidth',2)
xlabel('x')
ylabel('f_X')
legend('Sampled','Theory')
title('Asymmetric double-sided distribution')

%% Expected, mean and sample path
figure(2)
EX = (muS-0.5*sigma^2+lambda*(p/eta1-(1-p)/eta2))*t; % expected path
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
ylim([-1,1.2]);
title('Paths of a kou jump-diffusion process X = \mut + \sigmaW(t) + \Sigma_{i=1}^{N(t)} Z_i')

%% Variance = mean square deviation (Analytic expression for variance)

figure(3)
plot(t,(sigma^2+2*lambda*(p/eta1^2+(1-p)/eta2^2))*t,'b', t, var(X,0,2),'r',t,mean((X-transpose(EX)).^2,2),'g')
legend('Theory', "Sampled1", "Sampled2")
xlabel('t')
ylabel('Var(X) = E((X-E(X))^2)')
title('Kou jump-diffusion process X = \mut + \sigmaW(t) + \Sigma_{i=1}^{N(t)} Z_i: variance')

%% Probability density function at different times
figure(4)

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

print('-dpng','kjddensities.png')


