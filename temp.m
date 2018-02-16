% fun = @(x) exp(-x.^2).*log(x).^2;
% q = integral(fun,0,Inf)
L = 0.01;
q = 300E6
x = linspace(0,L,100);
k = 4;
T=100;
n =40;
fun = @(x) ((L^2*q*sin(pi.*x/L)/(pi^2*k))+T).*cos(2*n*pi.*x/(2*L));
top = integral(fun,0,L)

fun2 = @(x) (cos(n*pi.*x/(2*L))).^2;
bottom = integral(fun2,0,L)

bn = top/bottom