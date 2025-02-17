// Функция, задающая систему ОДУ:
function dx = f(t, x)
    dx = zeros(3,1);
    dx(1) = x(2);  
    dx(2) = x(3);  
    dx(3) = (10 - x(1) - 3.5*x(2) - 28*x(3)) / 8;
endfunction

// Метод Адамса-Мултона 1-го порядка (неявный метод) для системы ОДУ
function [T, X] = adams_moulton1(f, t0, tf, x0, h)
    N = int((tf - t0) / h);
    T = linspace(t0, tf, N+1);
    X = zeros(3, N+1);
    X(:,1) = x0;
    
    // Выполняем шаги метода
    for n = 1:N
        t_next = T(n+1);
        x_n = X(:,n);
        
        // Начальное приближение: шаг по методу Эйлера
        x_next = x_n + h * f(T(n), x_n);
        
        // Итерационное уточнение (5 итераций)
        for k = 1:5
            x_next = x_n + h * f(t_next, x_next);
        end
        
        X(:, n+1) = x_next;
    end
endfunction

// Заданные параметры задачи
t0 = 0;
tf = 100;
h = 0.1;  // Шаг по умолчанию



// Начальные условия: y(0) = 0, y'(0) = 0, y''(0) = 0
x0 = [0; 0; 0];

// Запуск численного метода
[T, X] = adams_moulton1(f, t0, tf, x0, h);

// Сохранение результатов в CSV-файл (столбцы: t и y(t)=x1)
data = [T' X(1,:)'];
csvWrite(data, "output.csv");

// Построение графика решения (plot2d, т.к. plot не работает в -nwni)
x = T;
y = X(1,:);

clf();  // Очистка графика перед построением
plot2d(x, y, style=2); // 2 — синий цвет

// Добавление подписей
xtitle("Решение ОДУ методом Адамса-Мултона 1-го порядка", "t", "y(t)");

// Сохранение графика в PNG через xs2png()
xs2png(gcf(), "output.png");

// Завершение работы скрипта
exit;
