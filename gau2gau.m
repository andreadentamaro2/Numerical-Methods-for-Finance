% Check numerically that Gaussians form a closed set with respect to
% the Fourier transform

N = 1024; % number of grid points (The algorithm performs better with a power of 2)
Dx = 0.1; % spacing of the x grid
x = Dx*(-N/2:N/2-1); % grid in x space
L = N*Dx; % truncation limit in x space
Dxi = 2*pi/L; % spacing of the xi grid; Nyquist relation: Dx*Dxi = 2*pi/N (specify the grid in Fourier space)
xi = Dxi*(-N/2:N/2-1); % grid in xi space
W = N*Dxi; % truncation limit in xi space; W = 2*pi/Dx

%%
a = 2;

% Analytical
fa = sqrt(a/pi)*exp(-a*x.^2); % Specify a Gaussian in x space
Fa = exp(-xi.^2/(4*a)); % Gaussian in xi space --> xi is the Fourier coniugated quantity of x 

% Numerical
Fn  = fftshift(ifft(ifftshift(fa)))*L; % Use ifft to do forward reverse fft
fn  = fftshift( fft(ifftshift(Fa)))/L; % Use fft to do down reverse fft
Fn1 = fftshift( fft(ifftshift(fa)))*Dx;
fn1 = fftshift(ifft(ifftshift(Fa)))/Dx;

close all
figure(1), clf, hold on
plot(x,real(fn),'r') % Real part of the inverse fuorier transform
plot(x,imag(fn),'g') % Imaginary part of the inverse fuorier transform
plot(x,fa,'y:') % Analytical part of the inverse fuorier transform
axis([-5 5 0 1])
xlabel('x')
ylabel('f')
legend('Re(fn)','Im(fn)','fa')

figure(2), clf, hold on
plot(xi,real(Fn),'b')
plot(xi,imag(Fn),'m')
plot(xi,Fa,'c:')
axis([-10 10 0 1])
xlabel('\xi')
ylabel('F')
legend('Re(Fn)','Im(Fn)','Fa')

%%
figure(3), clf, hold on
plot(x,real(fn1),'r')
plot(x,imag(fn1),'g')
plot(x,fa,'y:')
axis([-5 5 0 1])
xlabel('x')
ylabel('f')
legend('Re(fn1)','Im(fn1)','fa')

figure(4), clf, hold on
plot(xi,real(Fn1),'b')
plot(xi,imag(Fn1),'m')
plot(xi,Fa,'c:')
axis([-10 10 0 1])
xlabel('\xi')
ylabel('F')
legend('Re(Fn1)','Im(Fn1)','Fa')
