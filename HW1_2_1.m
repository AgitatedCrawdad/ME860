% clear 
% close
% Set maximum number of partial sums
nmax = 40;
c = 0.3;
k = 4.5;
Tw = 300;
Tw2 = 420;
rho = 10000;
alpha = k/(rho*c);
q = (300E6/(2/pi))*0.02;
L = 0.005;
% Initialize domain
x = linspace(-L,L,100);
% Plot initial condition

% Calculate temperature of rod at various times
t=[0.001 0.002 0.003 5];
count = 0;
for w=1:length(t)
T = 0;
for n=1:nmax
if n == 1
   Bn =(((4*(pi*k*L*Tw2-((L^3))*q))/((pi^2)*k)));
end
if n > 1 && mod(n,2)==0
    if mod(n,4) == 0
        Bn = -(16*(L^3)*q)/(((n^2)-1)*(pi^3)*k);
    else
        Bn = (16*(L^3)*q)/(((n^2)-1)*(pi^3)*k);
    end
else
    if mod(n+1,4) == 0
       Bn = -(4*L*Tw2)/(n*pi); 
    else
       Bn = (4*L*Tw2)/(n*pi);
    end
end

% Bn = 4*C/(n*pi);
% Bn = constant*((4*(pi*k*L*Tw-(2*(L^3))*q)/(n*(pi^2)*k)))/L;

% T =T +Bn*cos(n*pi*x/(2*L))*exp((-n^2*pi^2*alpha*t)/(4*L^2));
T =T +(Bn/L)*cos(n*pi*x/(2*L))*exp((-(n^2)*pi^2*alpha*t(w))/(4*L^2));

end

term2 = ((4*(L^2)*q*cos((pi*x)/(2*L))/((pi^2)*k)));
term3 = (-2*L*q*x/(pi*k));
term4 =  (2*(L^2)*q/(pi*k));
T = T + term2 + term3 + term4+Tw;
hold on
plot(x,T)
end
legend('t1','t2','t3','t=5'), xlabel('Length, x'), ylabel('Temperature')

