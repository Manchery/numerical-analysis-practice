a = 0; b = 20; h = 0.00001;
y0 = [0; -2];
n = floor((b-a)/h);

x = a + h*(0:n);
y_true = sol(x);

y_44 = myRungeKutta44(@f, y0, a, b, h);
y_24 = myRungeKutta24(@f, y0, a, b, h, 1e-12);

fprintf('max_ae y_44 %e\n', max_error(y_44, y_true));
fprintf('max_ae y_24 %e\n', max_error(y_24, y_true));

fprintf('mae y_44 %e\n', mae(y_44 - y_true));
fprintf('mae y_24 %e\n', mae(y_24 - y_true));

print_y(y_true, a, b, 4, 'y_true');
print_y(y_44, a, b, 4, 'y_44');
print_y(y_44-y_true, a, b, 4, 'y_44 - y_true');
print_y(y_24, a, b, 4, 'y_24');
print_y(y_24-y_true, a, b, 4, 'y_24 - y_true');

function res = f(x,y)
    u = y(1); v = y(2);
    du = -2000*u + 999.75*v + 1000.25;
    dv = u - v;
    res = [du; dv];
end

function y = sol(x)
    u = - 1.499875*exp(-0.5*x) + 0.499875*exp(-2000.5*x) + 1;
    v = - 2.99975*exp(-0.5*x) - 0.00025*exp(-2000.5*x) + 1;
    y = [u; v];
end

function print_y(y, a, b, n, label)
    disp(label);
    res = y(:, a+(b-a)/n);
    for i = 2:n
        res = [res, y(:, a+(b-a)*i/n)];
    end
    disp(res)
end