%% TEST 2021
% Candidate number: SRZZ4

%Question 1

%% Plot and sample the non-central chi-square distribution

clear all

% Define the parameters
d = 5; % degrees of freedom
lambda = 5; % non-centrality parameter
ngrid = 200;
nsample = 10^6;
a = 0;
b = 20;

% Define the grid
%x = linspace(0,20,100);
deltax = (b-a)/ngrid;
x = a:deltax:b;

%% Compute and plot the PDF and CDF
f = pdf('ncx2',x,d,lambda);
F = cdf('ncx2',x,d,lambda);

figure(1) % open a plot window
plot(x,f,'r',x,F,'b')
xlabel('x')
title('Non-central chi-square PDF and CDF with n=5 and d=2')
legend('PDF','CDF')
print('-dpdf','ncchisq.pdf')

%% (a) A normalised histogram that matches the PDF

X = ncx2rnd(d, lambda, 1, nsample);
figure(2)
x2 = [x-deltax/2 x(end)+deltax/2];
ho = histogram(X,x2,'normalization','pdf');
hold on
plot(x,f,'r','LineWidth',2)
xlim([a,b])
xlabel('x')
ylabel('f') 
legend('sampled','theory')
title('non-central chi squared pdf with d=5 and lambda=5')

%% (b) A scatter plot with 1000 points that cover uniformly the area under the PDF.

n = 1000 
X = ncx2rnd(d, lambda, 1, n);
U = rand(1,n);
f1 = pdf('ncx2',X,d,lambda);
Z = U.*f1;

figure(3);
plot(x,f, 'b', 'LineWidth',3)
hold on
scatter(X,Z,1,'filled')
xlim([a,b])
xlabel('x')
ylabel('U.*f1')
legend('scatter plot')
title('non-central chi-squared uniformly distributed in pdf')

