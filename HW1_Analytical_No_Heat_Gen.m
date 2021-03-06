% clear 
% close
% Set maximum number of partial sums
nmax = 40;
c = 300;
k = 4.5;
Tw = 300;
Tw2 = 420;
rho = 11000;
alpha = k/(rho*c);
q = (300E6/(2/pi))*0.02;
q1 = (300E6/(2/pi));
L = 0.01;
% Initialize domain
x = linspace(0,L,100);
% Plot initial condition

% Calculate temperature of rod at various times
% t=[1 2 3 4 5];
t = 0.001;
% t = linspace(0,5,10);
maxT=zeros(length(t),1);
count = 0;
for w=1:length(t)
T = 0;
for n=1:nmax
fun = @(x) ((Tw2-Tw)/2 + (((L^2)*q1*sin((pi.*x)/(L))/((pi^2)*k)))).*cos(n*pi.*x/(2*L));
top = integral(fun,0,L);

fun2 = @(x) (cos(n*pi.*x/(2*L))).^2;
bottom = integral(fun2,0,L);
Bn = top/bottom;

% Bn = 4*C/(n*pi);
% Bn = constant*((4*(pi*k*L*Tw-(2*(L^3))*q)/(n*(pi^2)*k)))/L;

% T =T +Bn*cos(n*pi*x/(2*L))*exp((-n^2*pi^2*alpha*t)/(4*L^2));
  T =T +Bn*cos(n*pi.*x/(2*L))*exp((-(n^2)*(pi^2)*alpha*t(w))/(4*L^2));



end

% term2 = (((L^2)*q*sin((pi*x)/(L))/((pi^2)*k)));
% term3 = (-L*q*x/(pi*k));
% term4 =  ((L^2)*q/(pi*k));
term2 = (((L^2)*q*sin(pi*x/L)/((pi^2)*k)));
term3 = (-L*q.*x/(pi*k));
term4 =  ((L^2)*q/(pi*k));
T = T + term2 + term3 + term4 + Tw;
% T = T + Tw;
hold on
plot(x,T,'--')
maxT(w) = max(T);
% plot(t(w),max(T),'--')
end
grid on
% plot(t,maxT,'--')
% xlabel('Time [s]')
xlabel('x [m]')
ylabel('Temperature [^{o}C]')

