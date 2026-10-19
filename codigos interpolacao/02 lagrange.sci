// ============================================================
// 02 - INTERPOLAÇÃO DE LAGRANGE
//      p(x) = soma_k y(k) * L_k(x)
//      L_k(x) = produto_{j<>k} (x - x(j)) / (x(k) - x(j))
// ============================================================

clear; clc;

// Pontos tabelados (x distintos) e ponto onde avaliar p(xp)
x  = [-1; 0; 2];
y  = [4; 1; -1];
xp = 1;

// Número de pontos (polinômio de grau n-1)
n = length(x);

// -------- Calcula cada L_k(xp) e acumula p(xp) --------
yp = 0;
for k = 1:n
    Lk = 1;
    for j = 1:n
        if j <> k
            Lk = Lk * (xp - x(j)) / (x(k) - x(j));
        end
    end
    L(k) = Lk;
    yp = yp + y(k) * Lk;
end

disp("=== INTERPOLAÇÃO DE LAGRANGE ===");
disp("Valores L_k(xp) na ordem k = 1..n:");
disp(L);
disp("p(xp) com xp =");
disp(xp);
disp(yp);
