// script2.sci

// Заданная матрица
A = [2 2 -3 3; 3 2 -1 1; 1 1 -2 2; 2 4 3 2];

// Функция LU-разложения с частичным выбором ведущего элемента
function [L, U, P] = LU_decomposition(A)
    [n, m] = size(A);
    if n ~= m then
        error("Матрица должна быть квадратной для LU-разложения.");
    end
    L = zeros(n, n);
    U = A;
    P = eye(n);
    for k = 1:n-1
        // Исправленный вызов max для Scilab
        [max_value, i_max] = max(abs(U(k:n, k)));
        i_max = i_max + k - 1;  // Индекс строки для перестановки
        if k ~= i_max then
            // Перестановка строк в U и P
            U([k, i_max], :) = U([i_max, k], :);
            P([k, i_max], :) = P([i_max, k], :);
        end
        for i = k+1:n
            L(i, k) = U(i, k) / U(k, k);
            U(i, k:n) = U(i, k:n) - L(i, k) * U(k, k:n);
        end
    end
    L = L + eye(n);  // Добавляем единичную матрицу на диагональ
endfunction


// Функция QR-разложения с использованием метода отражения
function [Q, R] = QR_decomposition(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;
    for k = 1:n
        x = R(k:m, k);
        e = zeros(length(x), 1);
        e(1) = 1;
        v = sign(x(1)) * norm(x) * e + x;
        v = v / norm(v);
        R(k:m, k:n) = R(k:m, k:n) - 2 * v * (v' * R(k:m, k:n));
        Q(k:m, :) = Q(k:m, :) - 2 * v * (v' * Q(k:m, :));
    end
endfunction

// Применение LU-разложения
[L_LU, U_LU, P_LU] = LU_decomposition(A);

// Применение QR-разложения
[Q_QR, R_QR] = QR_decomposition(A);

// Применение встроенных функций Scilab для LU и QR разложения
[L_builtin, U_builtin, P_builtin] = lu(A);
[Q_builtin, R_builtin] = qr(A);

// Запись результатов в CSV
csvWrite([L_LU; U_LU], "LU_result.csv");
csvWrite([Q_QR; R_QR], "QR_result.csv");

// Возврат к основным данным
exit;
