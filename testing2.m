% clear 
% close
% Set maximum number of partial sums
nmax = 100;
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
% t=[5 10 20 40 80];
% t = 5;
t = linspace(0,50,500);
maxT=zeros(length(t),1);
count = 0;
for n=1:nmax
fun = @(x) ((Tw2 + (((L^2)*q1*sin((pi.*x)/(L))/((pi^2)*k))))-Tw).*cos((2*n-1)*pi.*x/(2*L));
top = integral(fun,0,L);

fun2 = @(x) (cos((2*n-1)*pi.*x/(2*L))).^2;
bottom = integral(fun2,0,L);
Bn(n) = top/bottom;
end

for w=1:length(t)
T = 0;
for n=1:nmax
  T =T +Bn(n)*cos((2*n-1)*pi.*x/(2*L))*exp((-(4*n^2-4*n+1)*(pi^2)*alpha*t(w))/(4*L^2));
end

% for n=1:nmax
% fun = @(x) ((Tw2 + (((L^2)*q1*sin((pi.*x)/(L))/((pi^2)*k))))-Tw).*cos((2*n-1)*pi.*x/(2*L));
% top = integral(fun,0,L);
% 
% fun2 = @(x) (cos((2*n-1)*pi.*x/(2*L))).^2;
% bottom = integral(fun2,0,L);
% Bn = top/bottom;
% 
% % Bn = 4*C/(n*pi);
% % Bn = constant*((4*(pi*k*L*Tw-(2*(L^3))*q)/(n*(pi^2)*k)))/L;
% 
% % T =T +Bn*cos(n*pi*x/(2*L))*exp((-n^2*pi^2*alpha*t)/(4*L^2));
%   T =T +Bn*cos((2*n-1)*pi.*x/(2*L))*exp((-(4*n^2-4*n+1)*(pi^2)*alpha*t(w))/(4*L^2));
% end


% term2 = (((L^2)*q*sin((pi*x)/(L))/((pi^2)*k)));
% term3 = (-L*q*x/(pi*k));
% term4 =  ((L^2)*q/(pi*k));
% term2 = (((L^2)*q*sin(pi*x/L)/((pi^2)*k)));
% term3 = (-L*q.*x/(pi*k));
% term4 =  ((L^2)*q/(pi*k));
% T = T + term2 + term3 + term4 + Tw;
T = T + Tw;
% hold on

plot(x,T,'b--')
xlabel('x [m]')
ylabel('Temperature [^{o}C]')
legend(strcat('t = ', num2str(t(w))))
axis([0 0.01 300 1500])
grid on
MM(w) = getframe(gcf);
% maxT(w) = max(T);
% plot(t(w),max(T),'--')
% switch w
%     case 1
%         plot(x,T,'r--')
%     case 2
%         plot(x,T,'k--')
%     case 3
%         plot(x,T,'b--')
%     case 4
%         plot(x,T,'g--')
%     case 5
%         plot(x,T,'m--')
% end

end
movie(MM);
grid on
% plot(t,maxT,'b--')
% xlabel('Time [s]')
% USED TO MAKE A VIDEO
% myvideo = VideoWriter('myfile.avi');
% open(myvideo)
% writeVideo(myvideo, MM);
% close(myvideo)
