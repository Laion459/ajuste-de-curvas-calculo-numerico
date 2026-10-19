// ============================================================
// 03 - INTERPOLAÇÃO DE NEWTON (DIFERENÇAS DIVIDIDAS)
//      p(x) = d0 + d1*(x-x0) + d2*(x-x0)*(x-x1) + ...
//      d = diagonal da tabela de diferenças divididas
// ============================================================

clear; clc;

// Pontos tabelados (x distintos) e ponto onde avaliar p(xp)
x  = [-1; 0; 2];
y  = [4; 1; -1];
xp = 1;

// Número de pontos (polinômio de grau n-1)
n = length(x);

// -------- Monta a tabela de diferenças divididas --------
// Coluna 1 = ordem 0: D(i,1) = y(i)
for i = 1:n
    D(i, 1) = y(i);
end

// Coluna j = ordem j-1: D(i,j) = (D(i,j-1)-D(i-1,j-1)) / (x(i)-x(i-j+1))
for j = 2:n
    for i = j:n
        D(i, j) = (D(i, j-1) - D(i-1, j-1)) / (x(i) - x(i-j+1));
    end
end

// Coeficientes de Newton: diagonal da tabela (d0, d1, d2, ...)
for i = 1:n
    d(i) = D(i, i);
end

// -------- Avalia p(xp) na forma de Newton (Horner) --------
yp = d(n);
for k = n-1:-1:1
    yp = yp * (xp - x(k)) + d(k);
end

disp("=== INTERPOLAÇÃO DE NEWTON ===");
disp("Tabela de diferenças divididas D (coluna = ordem):");
disp(D);
disp("Coeficientes d = [d0; d1; d2; ...]:");
disp(d);
disp("p(xp) com xp =");
disp(xp);
disp(yp);
