nmax = 100;% Set maximum number of partial sums
k = 4.5; %Thermal Conductivity [W/mK]
rho = 11000; %Density of UO2 [kg/m^3]
c = 300; %Specific Heat Capcity of UO2  [J/kgK]
alpha = k/(rho*c); %Thermal Diffusivity [m^2/s]
q = (300E6/(2/pi))*0.02; %2% Peaking Power [W]
q1 = (300E6/(2/pi)); %Power at SS condition
Tw = 300; %Tw,b
Tw2 = 420; %Tw,a
% Initialize domain
L = 0.01;
x = linspace(0,L,100);
% Calculate temperature of rod at various times
% t=[5 10 20 40 80];
t = 5;
% t = linspace(0,5,11);
maxT=zeros(length(t),1);%Used to store max temps
for w=1:length(t)
T = 0;
for n=1:nmax
fun = @(x) ((Tw2 + (((L^2)*q1*sin((pi.*x)/(L))/((pi^2)*k))))-Tw).*cos((2*n-1)*pi.*x/(2*L)); %Top integral in eqn 20
top = integral(fun,0,L);

fun2 = @(x) (cos((2*n-1)*pi.*x/(2*L))).^2; %Bottom integral in eqn 20
bottom = integral(fun2,0,L); 
Cn = top/bottom;

  T =T +Cn*cos((2*n-1)*pi.*x/(2*L))*exp((-(4*n^2-4*n+1)*(pi^2)*alpha*t(w))/(4*L^2));%Eqn 19
end

T = T + Tw; %Adding initial wall temp back
hold on
plot(x,T,'b--')
maxT(w) = max(T);

end
grid on
% plot(t,maxT,'b--')
% xlabel('Time [s]')
xlabel('x [m]')
ylabel('Temperature [^{o}C]')

