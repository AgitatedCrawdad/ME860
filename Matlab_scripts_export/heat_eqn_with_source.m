clear close
% Set maximum number of partial sums
nmax = 40;
g = 1;
alpha=1;
C=1;
% Initialize domain
x = linspace(0,3,1000);
% Plot initial condition
for i = 1:length(x)
T (i)= 0;
%end 
end
plot(x,T) 

hold on
% Calculate temperature of rod at various times
t=[0.2,1,3,100];
for k = 1:4
T = -g*x.^2/2 + g*3*x/2;
for n=1:2:nmax
%Bn = 4*C/(n*pi);
Bn = -4*3^2*g/(n^3*pi^3);
T =T+ Bn*sin(n*pi*x/3)*exp(- n^2*pi^2*alpha*t(k)/9);
end
plot(x,T)
end
legend('t=0','t=0.2','t=1','t=3','t=100'), xlabel('Length along rod, x'), ylabel('Temperature, C')

