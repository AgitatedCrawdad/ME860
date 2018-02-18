% clear 
% close
% Set maximum number of partial sums
nmax = 40;
C=300;
alpha=1;
L = 0.01;
% Initialize domain
x = linspace(0,L,100);
% Plot initial condition
Tw = 420;
Tinf = 300;
k = 4.25;
h = 25000;
q = 300E6/(2/pi);
% T = Tw + ((4*(L^2)*q*cos((pi*x)/(2*L))/((pi^2)*k))); %From -L to L
T = Tw + (((L^2)*q*sin((pi*x)/(L))/((pi^2)*k)));       %From 0 to L

% T =Tw-((4*(L^2)*q*cos((pi*x)/(2*L))/((pi^2)*k)))-(h*(x+L)*(Tw-Tinf)/k)-((2*L*q*(x+L)/(pi*k)));
plot(x,T)
grid on
xlabel('x [m]')
ylabel('Temperature [^{o}C]')