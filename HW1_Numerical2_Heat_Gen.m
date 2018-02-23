% clear 
% close all
% Initialize number of nodes and constants
tic
N = 100;

L = 0.01;
rho = 11000;
k = 4.5;
c = 300;
alpha = k/(rho*c);
q = 300E6/(2/pi);
h = L/(N-1);
dt = 0.0001;
% Initialize domain
x = linspace(0,L,N);
x1 = x;
% Initialize and plot heat distribTtion at t = 0
T = ones(1,N)*500; 
T = T;
% hold on
% Calculate heat at specified times
% t = [1,1.1,6,8];
t = 5;
jmax = t/dt;
k = 4.5;
T(1) = 420;
T(length(T)) = 420; 
error = 1;
count = 0;
while error > 0.00001
Told = T;
for i=2:N-1    
    T(i) = (Told(i-1) + Told(i+1))/2 + (h^2)*q*sin((pi*i*h)/L)/(2*k);
% T(i) = Told(i) + alpha*dt((Told(i-1)- 2*Toldi*(i)+Told(i+1))/(h^2) + q*sin((pi*i*h)/L));
end
error = max(abs(T-Told));
end
% plot(x,T)
grid on
% xlabel('x')
% ylabel('T')


% clear 
% close all
% Initialize number of nodes and constants

L = 0.01;
rho = 11000;
k = 4.5;
c = 300;
alpha = k/(rho*c);
q = (300E6/(2/pi))*0.02;
h = L/(N-1);
dt = 0.0037;
% Initialize domain
x = linspace(0,L,N);
x1 = x;
% Initialize and plot heat distribTtion at t = 0
% hold on
% Calculate heat at specified times
% t = [5,10,20,40,80];
% t = 5;
t = linspace(0,5,11);
jmax = t/dt;
maxT=zeros(length(t),1);
w = 1;
T(length(T)) = 300;
for j = 1:jmax(length(t)) 
    Told = T;
for i=2:N-1    
%     T(i) = Told(i) + alpha*dt*((Told(i-1)-2*Told(i)+Told(i+1))/(h^2) + q);
%       T(i) = Told(i) + alpha*dt*((Told(i-1)-2*Told(i)+Told(i+1))/(h^2) + q*sin((pi*i*h)/L)/k);
      T(i) = Told(i) + alpha*dt*((Told(i-1)-2*Told(i)+Told(i+1))/(h^2));
end
T(1)=T(2);
if j == jmax(w)||jmax(w)==0 
maxT(w) = max(T);
% plot(x,T)
% 
% switch w
%     case 1
%         plot(x,T,'r')
%     case 2
%         plot(x,T,'k')
%     case 3
%         plot(x,T,'b')
%     case 4
%         plot(x,T,'g')
%     case 5
%         plot(x,T,'m')
% end

% xlabel('x')
% ylabel('T')
w = w+1;
end
plot(x,T,'r')
MM(j) = getframe; 
end
movie(MM)
toc

% plot(t,maxT,'r')
% legend('t1')
% xlabel('Length, x')
% ylabel('Temperature')