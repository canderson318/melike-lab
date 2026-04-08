clear; clc;

%%
f = @(t,y) y.^2; 

Y0 = -3; % y at t=0 
[t,y] = ode45(f, [0 5] , Y0 ); 
plot(t,y, 'o')
hold on 
syms y(t)
eq = diff(y) == y^2;
sol = dsolve(eq, y(0) == -3);
fplot(sol, [0 5])
