%% {Uniform function or Rectangual pulse} to sinc function

N = 1024; % number of grid points
Dx = 0.01; % spacing of the x grid 
x = Dx*(-N/2:N/2-1); % grid in x space
Lx = N*Dx; % truncation limit in x space
Dxi = 2*pi/Lx; % spacing of the xi grid; Nyquist relation: Dx*Dxi = 2*pi/N
xi = Dxi*(-N/2:N/2-1); % grid in xi space
Lxi = N*Dxi; % truncation limit in space; W = 2*pi/Dx

%%
a = 2;
fa = (x>=-a).*(x<=a);
Fa = 2*sin(a*xi)./xi;
Fa(end/2+1) = 2*a

Fn  = fftshift(ifft(ifftshift(fa)))*Lx;
fn  = fftshift( fft(ifftshift(Fa)))/Lx;
Fn1 = fftshift( fft(ifftshift(fa)))*Dx;
fn1 = fftshift(ifft(ifftshift(Fa)))/Dx;

close all

figure(1), clf, hold on
plot(x,real(fn),'r')
plot(x,imag(fn),'g')
plot(x,fa,'b')
axis([-5 5 -0.5 1.5])
xlabel('x')
ylabel('f')
legend('Re(fn)','Im(fn)','fa')

figure(2), clf, hold on
plot(xi,real(Fn),'b')
plot(xi,imag(Fn),'m')
plot(xi,Fa,'b:')
axis([-100 100 -1 2*a])
xlabel('\xi')
ylabel('F')
legend('Re(Fn)','Im(Fn)','Fa')

figure(3), clf, hold on
plot(x,real(fn1),'r')
plot(x,imag(fn1),'g')
plot(x,fa,'b')
axis([-5 5 -0.5 1.5])
xlabel('x')
ylabel('f')
legend('Re(fn1)','Im(fn1)','fa')

figure(4), clf, hold on
plot(xi,real(Fn1),'b')
plot(xi,imag(Fn1),'m')
plot(xi,Fa,'b:')
axis([-100 100 -1 2*a])
xlabel('\xi')
ylabel('F')
legend('Re(Fn1)','Im(Fn1)','Fa')