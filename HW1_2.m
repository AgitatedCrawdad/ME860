clear 
close
% Set maximum number of partial sums
nmax = 40;
C=300;
alpha=1;
L = 0.005;
% Initialize domain
x = linspace(-L,L,100);
% Plot initial condition
Tw = 300;
Tinf = 300;
k = 4.5;
h = 25000;
q = 300E6/(2/pi);
term1 = ((4*(L^2)*q*cos((pi*x)/(2*L))/((pi^2)*k)));
term2 = (-2*L*q*x/(pi*k));
term3 = (2*(L^2)*q/(pi*k));
T = Tw + term1 + term2 + term3;
% T =Tw-((4*(L^2)*q*cos((pi*x)/(2*L))/((pi^2)*k)))-(h*(x+L)*(Tw-Tinf)/k)-((2*L*q*(x+L)/(pi*k)));
plot(x,T)
grid on