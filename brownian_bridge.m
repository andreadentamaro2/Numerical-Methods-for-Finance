%% Monte Carlo simulation of the Brownian bridge
% dX = (b-X)/(T-t)*dt + sigma*dW

% Define parameters and time grid
clear all % clear all variables from memory
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
sigma = 1; % volatility
a = 0.8; % initial value
b = 1; % final value

%% Monte Carlo method 1
% Update noice and drift part with each update 
%Allocate and initialise all paths
X = [a*ones(1,npaths); zeros(nsteps-1,npaths); b*ones(1,npaths)];

% Compute the Brownian bridge with Euler-Maruyama
for i = 1:nsteps-1
     X(i+1,:) = X(i,:) + (b-X(i,:))/(nsteps-i+1) + sigma*randn(1,npaths)*sqrt(dt);
end

%% Monte Carlo method 2

%Compute the increments of driftless arithmetic Brownian motion
%dW = sigma*sqrt(dt)*randn(nsteps,npaths);

% Accumulate the increments of arithmetic Brownian motion
%W = cumsum([a*ones(1,npaths); dW]); % a could be unequal 0

% Compute the Brownian bridge with X(t) = W(t) + (b-W(T))/T*t
%X = W + repmat(b-W(end,:),nsteps+1,1)/T.*repmat(t',1,npaths);

%% Expected, mean and sample paths
figure(1)
close all
EX = a + (b-a)/T*t; % expected path
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
%sdevmax = sigma*sqrt(T)/2;
%ylim([(a+b)/2-4*sdevmax (a+b)/2+4*sdevmax])
title('Brownian bridge dX = ((b-X)/(T-t))dt + \sigmadW')
print('-dpng','bbpaths.png')

%% Variance = mean square deviation (Analytic expression for variance)

figure(2)
plot(t, (T-t)/T.*t,'b', t, var(X,0,2), 'r', t,mean((X-transpose(EX)).^2,2),'g');
legend('Theory','Sample 1', 'Sample 2')
xlabel('t')
ylabel('Var(X) = E((X-E(X))^2)')
title('Brownian Bridge: variance')

%% Plot the probability density function at different times
figure(3)

subplot(3,1,1)
histogram(X(20,:),0:0.035:3.5,'normalization','pdf');
ylabel('f_X(x,0.15)')
xlim([0,1.5])
ylim([0,7])
title('Brownian bridge: PDF at different times')

subplot(3,1,2)
histogram(X(80,:),0:0.035:3.5,'normalization','pdf');
xlim([0,1.5])
ylim([0,7])
ylabel('f_X(x,0.4)')

subplot(3,1,3)
histogram(X(end-1,:),0:0.035:3.5,'normalization','pdf');
xlim([0,1.5])
ylim([0,7])
xlabel('x')
ylabel('f_X(x,1)')
print('-dpng','gbpdensities.png')

figure(4)
s = 4
T = 10
nSteps = 100

tRange = linspace(s, 10, nSteps)

[z, t] = meshgrid(tRange, tRange)

%%

covFunc = @(zArr, tArr) ((min(zArr, tArr) - s) * (T - max(zArr, tArr))) / (T - s)

covAll = covFunc(z, t)


surf(t,z,covAll)


