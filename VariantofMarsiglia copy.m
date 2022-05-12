%% Variant of Marsiglia
 clear all

% Define the parameter
mu = 0.2;
sigma = 0.1;    
n = 10000

w = zeros(n);
for i = 1:n
    w = 1;
    while w >= 1
        U1 = rand();
        U2 = rand();
        V1 = (2*U1)-1;
        V2 = (2*U2)-1;
        w = V1^2 + V2^2;
    end
    N(i) = V2*sqrt((-2*log(w))/w);
end



%% Plotting the analytical probability density function
a = -0.2; %left truncation
b = 0.6; %right truncation
x = linspace(a,b,200); %grid with 200 equal steps
nbins = 100

%analytical normal PDF with x as an argument
f1 = 1/(sqrt(2*pi)*sigma)*exp(-((x-mu)/sigma).^2/2);

%analytical normal PDF with Z as an argument
f2 = 1/(sqrt(2*pi)*sigma)*exp(-((N-mu)/sigma).^2/2);

%New uniform random values
U=rand(1,n);

figure(1)
plot(x,f1); %plotting analytical PDF
hold on
plot(N,f2.*U,'.'); %covering area underneeth with random points
xlim([a b])
xlabel('x') % label of the x axis
legend('PDF', 'Random Points')
title('Normal distribution with \mu = 0.2 and \sigma = 0.1')


figure(2)
histogram(N,'Normalization','pdf')