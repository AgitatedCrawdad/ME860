% clear 
% close
% Set maximum number of partial sums
nmax = 40;
C=1;
alpha=1;
L = 3
% Initialize domain
x = linspace(0,3,1000);
% Plot initial condition
for i = 1:length(x)
T(i) = C;
end
plot(x,T) 
hold on
% Calculate temperature of rod at various times
t=[0.2,1.5,4,10];
for k = 1:4
T = 0;
for n=1:2:nmax
Bn = 4*C/(n*pi);
T =T+ Bn*sin(n*pi*x/L)*exp(- n^2*pi^2*alpha*t(k)/L^2);
end
plot(x,T)
end
legend('t=0','t=0.2','t=1.5','t=4','t=10'), xlabel('Length, x'), ylabel('Temperature')

