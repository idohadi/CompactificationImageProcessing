l = 20; 
m = -10; 


N = 10^3; 
epsi = 10^-3; 

x = linspace(-1+epsi, 1-epsi, 2*N); 
t = tic();
y = evalAlegendre(l, m, x);
t = toc(t);
disp(['Runtime (sec): ', num2str(t)]);
plot(x, y)

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
