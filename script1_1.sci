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

// Сетка времени для встроенной функции ode
T_grid = linspace(t0, tf, int((tf - t0) / h) + 1);

// Решение с помощью встроенной функции ode
[X_builtin, T_builtin] = ode(x0, t0, T_grid, f);

// Преобразование в формат (N, 3) для корректного объединения данных
X_builtin_transposed = X_builtin';

// Проверка размерностей и запись данных
if size(T_rk4, "c") == size(T_builtin, "c") then
    data = [T_rk4' X_rk4(1,:)' X_builtin_transposed(:,1)];
    csvWrite(data, "output.csv");
else
    printf("Размерности временных шагов не совпадают: Рунге-Кутта — %d, ode — %d\n", size(T_rk4, "c"), size(T_builtin, "c"));
end

exit;
