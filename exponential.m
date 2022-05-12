% Plot and sample the exponential distribution

% Define the parameters
mu = 2/3;
nsample = 10^6;

% Define the grid
%x = linspace(0,5,101);
deltax = 0.05
x = 0:deltax:5;

% Compute the PDF and the CDF
% Unlike the normal distribution we can find both PDF and CDF analytically
f = 1/mu*exp(-x/mu);
F = 1-exp(-x/mu);
%f = pdf('Exponential',x,mu);
%F = cdf('Exponential',x,mu);

% Open a plot window
close all
figure(1)
plot(x,f,'r',x,F,'b');

% Label the x axis
xlabel('x')

% Insert the legend
legend('PDF','CDF')

% Set the title
title('Exponential distribution with \mu = 2/3')

% Print the figure to a file
print('-dpdf','exponential.pdf')

% Sample an exponential distribution
U = rand(nsample,1); % standard uniform random numbers
X = -mu*log(U); % transformation to exponential random numbers
% The transformation is done by hand because it is easy to find the inverse
% of the CDF, which is X = -mu*log(U)

% Bin the random variables in a histogram and normalise it
dx = x(2)-x(1); % bin width
[h,xx] = hist(X,x+dx/2);
h = h/(nsample*dx); % normalisation

figure(2)
plot(xx,h,'b',x,f,'r')
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('Exponential distribution with \mu = 0.2')

figure(3)
bar(xx,h)
hold on
plot(x,f,'r')
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('Exponential distribution with \mu = 0.2')

figure(5)
x2 = [x+dx/2 x(end)+dx/2]; % bin edges
ho = histogram(X,x2,'normalization','pdf'); % histogram object (MOST IMPORTANT COMMAND)
hold on
plot(x,f,'r')
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('Exponential distribution with \mu = 0.2')
 