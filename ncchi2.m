%% TEST 1 2021

%% Plot the non-centra chi-squared distribution

% Assign Parameters
n = 5; % degrees of freedom
d = 5; % non-centrality parameter 

% Assign x range
x = linspace(0,20,100);

%% Compute and plot the PDF and CDF
fncchi2 = pdf("ncx2", x, n, d);
Fncchi2 = cdf("ncx2", x, n, d);

% Make the plot 
figure(1)
plot(x, fncchi2, "r", x, Fncchi2, "b")
xlabel("x")
title("Non central Chi-Square pdf and cdf when n=5 and d=2")
legend("Non Central Chi-Square pdf, Non Central Chi-Square cdf")

%% Sample the non-central chi-squaree distribution

%(a.)

d = 5;
lambda = 5;
nsample =20000;
ngrid = 100;
a = 0; % left truncation
b = 20; % right truncation
deltax = (b-a)/ngrid;


%X = icdf("ncx2", rand(1,nsample), d, lambda); % 50 times slower
X = ncx2rnd(d, lambda, 1, nsample);

figure(2)
x2 = [x-deltax/2 x(end)+deltax/2]; % bin edges 
ho = histogram(X, x2, "Normalization","pdf"); % histogram object
hold on
plot(x,fncchi2,"r", "Linewidth", 2);
xlim([a b]);
xlabel("x")
ylabel("f")
legend("Sampled, Theory")
title("Non-central chi-square PDF with d=5 and \lambda=5")

figure(3)
h2 = ho.Values;
plot(x, h2, "b", x, fncchi2, "r--", "LineWidth", 2)
xlim([a,b])
xlabel("x")
ylabel("f")
legend("Sampled", "Theory")
title("Non-central chi-square PDF with d=5 and \lambda=5")


% (b)

%New uniform random values
U=rand(1,n);

figure(4)
plot(x, fncchi2);
hold on




