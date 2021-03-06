% clear 
% close all
% Initialize number of nodes and constants
N = 100;

L = 0.01;
rho = 11000;
k = 4.5;
c = 300;
alpha = k/(rho*c);
q = 300E6/(2/pi)*0.02;
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
t = [1,1.1,6,8];
jmax = t/dt;
k = 4.5;
T(1) = 300;
T(length(T)) = 300; 
error = 1;
count = 0;
while error > 0.00001
Told = T;
for i=2:N-1    
%     T(i) = (Told(i-1) + Told(i+1))/2 + (h^2)*q*sin((pi*i*h)/L)/(2*k);
T(i) = (Told(i-1) + Told(i+1))/2 + (h^2)*q*sin((pi*i*h)/L)/(2*k);
% T(i) = Told(i) + alpha*dt((Told(i-1)- 2*Toldi*(i)+Told(i+1))/(h^2) + q*sin((pi*i*h)/L));
end
T(1) = T(2);
error = max(abs(T-Told));
end
plot(x,T)
grid on
xlabel('x')
ylabel('T')
