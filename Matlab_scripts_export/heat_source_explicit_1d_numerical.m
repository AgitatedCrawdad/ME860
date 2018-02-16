clear 
close all
% Initialize number of nodes and constants
N = 10;
alpha=1;
L = 3;

g = 1;
h = L/(N-1); dt = 0.0001;
% Initialize domain
x = linspace(0,L,N);
x1 = x;
% Initialize and plot heat distribTtion at t = 0
T = zeros(1,N); 
hold on
% Calculate heat at specified times
t = [0.02,0.05,0.1,1];
jmax = t/dt;
k = 1;
for j = 1:jmax(4) 
    Told = T;
for i=2:N-1    
T(i) = Told(i) + alpha*dt/h^2*(Told(i-1)- 2*Told(i)+Told(i+1)) + g*dt;
end
% clf
% plot(x,T)
% xlim([0 3])
% ylim([0 2])
% pause(0.00000001)

if j == jmax(k)
plot(x,T)
xlabel('x')
ylabel('T')
k = k+1;
end
end