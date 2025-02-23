// script1_1.sci

// Функция, задающая систему ОДУ:
function dx = f(t, x)
    dx = zeros(3, 1);
    dx(1) = x(2);
    dx(2) = x(3);
    dx(3) = (10 - x(1) - 3.5*x(2) - 28*x(3)) / 8;
endfunction

// Метод Рунге–Кутта 4-го порядка
function [T, X] = runge_kutta4(f, t0, tf, x0, h)
    N = int((tf - t0) / h);
    T = linspace(t0, tf, N+1);
    X = zeros(length(x0), N+1);
    X(:,1) = x0;
    for n = 1:N
        t_n = T(n);
        x_n = X(:, n);
        k1 = f(t_n, x_n);
        k2 = f(t_n + h/2, x_n + (h/2)*k1);
        k3 = f(t_n + h/2, x_n + (h/2)*k2);
        k4 = f(t_n + h, x_n + h*k3);
        X(:, n+1) = x_n + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
endfunction

// Параметры
t0 = 0;
tf = 100;

try
    h_str = getenv("STEP");
catch
    h_str = "";
end

if h_str == "" then
    h = 0.1;
else
    h = evstr(h_str); // Преобразуем строку в число
end

// Начальные условия
x0 = [0; 0; 0];

// Решение ОДУ методом Рунге–Кутта 4-го порядка
[T_rk4, X_rk4] = runge_kutta4(f, t0, tf, x0, h);

// Сетка времени
N = int((tf - t0) / h);
T_grid = linspace(t0, tf, N+1);

// Решение с помощью встроенной функции ode (аналог ode45)
X_builtin = ode(x0, t0, T_grid, f);

// Сохранение данных в CSV
data = [T_rk4' X_rk4(1,:)' X_builtin(:,1)];
csvWrite(data, "output.csv");

exit;
