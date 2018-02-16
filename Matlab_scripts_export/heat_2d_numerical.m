clear 
close
% Initialize number of nodes
M = 10;
N = 10;
% Initialize function
u = zeros(M,N);
% Initialize boundary conditions 
u(M,:) = 100;
% Calculate temperature distribution
flag = 1;
while flag == 1
flag = 0; 
uold = u;
for i = 2:M-1
for j = 2:N-1
u(i,j) = (u(i-1,j) + uold(i+1,j) + u(i,j-1) + uold(i,j+1))/4;
error = abs((u(i,j)-uold(i,j))/uold(i,j));
if error > 0.001
flag = 1;
end
end
end
end
% Plot temperature distribution
[x,y] = meshgrid(0:1/(M-1):1 , 0:1/(N-1):1);
surface(x,y,u)
colorbar
% colormap gray
title('2D steady state, Gauss-Seidel, 40 nodes')
xlabel('x')
ylabel('y')