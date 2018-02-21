L = 0.01; %Length of the domain [m]
% Initialize domain
x = linspace(0,L,100);
% Plot initial condition
Tw = 420; %Temperature at the wall [Celsius]
k = 4.5; %Thermal Conductivity [W/mK]
q = 300E6/(2/pi); %Converts Average Heat Gen. To Peak Heat Gen
T = Tw + (((L^2)*q*sin((pi*x)/(L))/((pi^2)*k)));       %From 0 to L
plot(x,T)
grid on
xlabel('x [m]')
ylabel('Temperature [^{o}C]')