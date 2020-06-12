l = 200; 
m = -50; 


N = 10^3; 
epsi = 10^-3; 

x = linspace(-1+epsi, 1-epsi, 2*N); 
y = evalAlegendre(l, m, x);
plot(x, y)

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
