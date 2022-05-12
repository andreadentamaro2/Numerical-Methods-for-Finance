%% Accept-reject method for normal from Laplace random variables

clear all
close all

nsample =1000000;
c = (2/sqrt(2*pi))*exp(0.5); % max(f/g) set x=1
dx = 0.1; % grid spacing
xmax = 4; % grid bound
x = -xmax:dx:xmax; % grid

% Sample the standard Laplace or double-sided exponential distribution
U1 = rand(1,nsample);
L = log(2*U1).*(U1<0.5)-log(2*(1-U1)).*(U1>=0.5); %Quantile function

% Sample the normal distribution using the accept-reject algorithm
g = 0.5*exp(-abs(L));
f = 1/(sqrt(2*pi))*exp(-0.5*L.^2);
U2 = rand(1,nsample);
N = L(U2*c.*g<=f);
U = U2(U2*c.*g<=f);

% Output to console
length(N)/nsample % acceptance ratio
format long
c % analytical value
max(f./g) % numerical check

figure(1) %Sampled normal distribution function
x2 = [x-dx/2 x(end)+dx/2]; % bin edges
histogram(N,x2,'normalization','pdf');
hold on
fx = 1/sqrt(2*pi)*exp(-(x.^2)/2);
gx = 0.5*exp(-abs(x));
plot(x,fx,x,c*gx,'g') % We multiply c time g otherwise the parabola would touch the histogram 
xlabel('x')
legend('Sampled f(x)','Theoretical f(x)','Majorant function cg(x)')
title('Standard normal distribution using the accept-reject algorithm')

figure(2)
plot(x,x.^2-2*x+1,x,fx./gx,x,c*ones(size(x)),'--g')
xlim([0 3])
xlabel('x')
legend('x^2-2x+1','f/g','c = (2e/pi)^{1/2}','Location','northwest')
title('x^2-2x+1 = 0 where f/g = max = c')


