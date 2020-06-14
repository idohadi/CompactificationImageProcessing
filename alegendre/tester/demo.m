l = 10; 
m = -5; 


N = 5*10^3+1; 

x = linspace(-1, 1, N); 
t = tic();
y = evalAlegendre(l, m, x);
t = toc(t);
disp(['Runtime (sec): ', num2str(t)]);
plot(x, y)

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
