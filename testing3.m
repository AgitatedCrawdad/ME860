%%%%%USE T FROM PART A AS INITIAL CONDITION%%%%%
k = 4.5; %Thermal Conductivity [W/mK]
rho = 11000; %Density of UO2 [kg/m^3]
c = 300; %Specific Heat Capcity of UO2  [J/kgK]
alpha = k/(rho*c); %Thermal Diffusivity [m^2/s]
q = (300E6/(2/pi))*0.02; %2% Peaking Power [W]
h = L/(N-1); %Delta x
dt = 0.001; %Delta t
% Initialize domain
x = linspace(0,L,N);
% Calculate heat at specified times
% t = [5,10,20,40,80]; %Used to calculate at multiple times
t = 20; %Used to calculate at one time
% t = linspace(0,5,11);%Used to calculate at multiple times
jmax = t/dt; %Used to find number of time steps
maxT=zeros(length(t),1);%Used to max temp.
w = 1;
T(length(T)) = 300; %Boundary Condition of Wall Temp = 300C
for j = 1:jmax(length(t)) 
    Told = T;
for i=2:N-1    
      T(i) = Told(i) + alpha*dt*((Told(i-1)-2*Told(i)+Told(i+1))/(h^2) + q*sin((pi*i*h)/L)/k); %FDE with heat generation 
%       T(i) = Told(i) + alpha*dt*((Told(i-1)-2*Told(i)+Told(i+1))/(h^2)); %FDE without heat generation 
end
T(1)=T(2); %Boundary Condition used to model insulation at x=0
if j == jmax(w)||jmax(w)==0 
maxT(w) = max(T);
w = w+1;
% plot(x,T,'r')
end
end
plot(x,T,'r')
