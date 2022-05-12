%% Fibonacci random number generator (Knuth)
% Seydel, Course Notes, Chapter 2, pages 211-212

% Parameters
mu = 5;
nu = 17;
M = 714025;
a = 1366;
b = 150889;


seed = 12345;

nsample = 10000;
nbins = 100

% Linear congruential generator
U = zeros(1,nu);
U(1) = mod(seed,M);
for i = 2:nsample
    U(i) = mod(a*U(i-1)+b,M);
end

% Fibonacci Generator
for i = nu+1:nsample
    U(i) = mod(U(i-nu)-U(i-mu),M);
end

U = U/M

%% Plot

% 2D scatter plot
figure(1)
plot(U(1:nsample-1),U(2:nsample),".")
xlabel('U_{i-1}')
ylabel('U_i')
title('Scatter plot')

% Probability density function
figure(2)
[h,x] = hist(U,100);
h = h*nbins/nsample
plot(x,h,x,ones(size(x)))
xlim([-0.2 1.2])
ylim([0 2])
xlabel("x")
ylabel("f")
legend("Sampled","Theory")
title("Uniform Distribution")


   

