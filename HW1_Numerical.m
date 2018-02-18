% clear 
% close all
% Initialize number of nodes and constants
N = 100;

L = 0.01;
rho = 11000;
k = 4.5;
c = 300;
alpha = k/(rho*c);
q = (300E6/(2/pi));
h = L/(N-1);
dt = 0.0001;
% Initialize domain
x = linspace(0,L,N);
x1 = x;
% Initialize and plot heat distribTtion at t = 0
T = ones(1,N)*500; 
T = T;
hold on
% Calculate heat at specified times
t = [1,1.1,6,10];
jmax = t/dt;
w = 1;
T(1) = 420;
T(length(T)) = 420;
for j = 1:jmax(4) 
    Told = T;
for i=2:N-1    
    T(i) = Told(i) + alpha*dt*((Told(i-1)-2*Told(i)+Told(i+1))/(h^2) + q*sin((pi*i*h)/L)/k);
% T(i) = Told(i) + alpha*dt((Told(i-1)- 2*Toldi*(i)+Told(i+1))/(h^2) + q*sin((pi*i*h)/L));
end
% clf
% plot(x,T)
% xlim([0 3])
% ylim([0 2])
% pause(0.00000001)

if j == jmax(w)
plot(x,T)
xlabel('x')
ylabel('T')
w = w+1;
end
end