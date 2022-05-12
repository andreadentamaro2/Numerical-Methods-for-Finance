%% Correlated Random Numbers

clear all

% Define parameters
n = 1000;
mu1 = 0.2; mu2 = 0.25;
sigma1 = 0.1;   sigma2 = 0.15;

% Generate normal random vectors
X1 = randn(1,n);
X2 = randn(1,n);
cov = mean((X1-mu1).*(X2-mu2));

% Expected mean
EX1 = mu1+(X1-mu2);
EX2 = mu2+(X2-mu2);

rho = cov/(sigma1*sigma2);

%% Correlated Random Numbers

clear all

% Define the parameters
sigma = 0.5;
n = 10000

X = sigma*randn(1,n);
Y = sigma*randn(1,n);

rho = 0.2;
A = X;
B = sqrt(rho^2)*X+sqrt(1-rho^2)*Y;
[std(A), std(B)]
corrcoef(A, B)







