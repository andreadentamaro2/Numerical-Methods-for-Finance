%% TEST 2021
% Candidate number: SRZZ4

% Question 2

%% Compute the Black-Scholes price of a Lookback Option

clear all

% Market parameters
npaths = 20000;
T = 0.5; % maturity
S0 = 1; % spot price
K = 1; % strike price
r = 0.1; % risk-free interest rate
q = 0; % dividend rate
S0 = 1; % initial price
nsteps = 50 ;
dt = T/nsteps;

% Model parameter
sigma = 0.3; % volatility

% Risk-neutral measure
muRN = r-q-0.5*sigma^2; % drift coefficient of the arithmetic Brownian motion

%% Monte Carlo
%(a.) The call price of a fixed-strike lookback option on the maximum
tic
% Arithmetic Brownian motion X(T) = log(S(T)/S(0)) at time T
dX = muRN*dt + sigma*randn(npaths,nsteps)*sqrt(dt);

% Accumulate the increments
X = [zeros(npaths,1) cumsum(dX,2)];

% Transform to geometric Brownian motion
S = S0*exp(X);

max_path = max(S,[],2);
VcMCb = exp(-r*T)*mean(max(max_path-K,0));
VcMCb
toc

% (b.) The put price of a fixed-strike lookback option on the minimum
tic
% Arithmetic Brownian motion X(T) = log(S(T)/S(0)) at time T
dX = muRN*dt + sigma*randn(npaths,nsteps)*sqrt(dt);

% Accumulate the increments
X = [zeros(npaths,1) cumsum(dX,2)];

% Transform to geometric Brownian motion
S = S0*exp(X);

min_path = min(S,[],2);
VpMCb = exp(-r*T)*mean(max(K-min_path,0));
VpMCb
toc




