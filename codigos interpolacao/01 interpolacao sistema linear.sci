// ============================================================
// 01 - INTERPOLAÇÃO VIA SISTEMA LINEAR (VANDERMONDE)
//      Monta V * a = y e resolve por Eliminação de Gauss
//      p(x) = a0 + a1*x + a2*x^2 + ... (caso base: linear/quadrática)
// ============================================================

clear; clc;

// Pontos tabelados (x distintos) e ponto onde avaliar p(xp)
x  = [-1; 0; 2];
y  = [4; 1; -1];
xp = 1;

// Número de pontos (polinômio de grau n-1)
n = length(x);

// -------- Monta a matriz de Vandermonde V --------
// Linha i: [1, x(i), x(i)^2, ..., x(i)^(n-1)]
for i = 1:n
    for j = 1:n
        V(i, j) = x(i)^(j-1);
    end
end

// Matriz aumentada [V | y]
Ab = [V, y];

// -------- FASE 1: Eliminação para frente (igual ao 01 gauss) --------
for k = 1:n-1
    for i = k+1:n
        fator = Ab(i, k) / Ab(k, k);
        Ab(i, k:n+1) = Ab(i, k:n+1) - fator * Ab(k, k:n+1);
    end
end

// -------- FASE 2: Substituição retroativa --------
a(n) = Ab(n, n+1) / Ab(n, n);
for i = n-1:-1:1
    soma = 0;
    for j = i+1:n
        soma = soma + Ab(i, j) * a(j);
    end
    a(i) = (Ab(i, n+1) - soma) / Ab(i, i);
end

// -------- Avalia p(xp) = a0 + a1*xp + a2*xp^2 + ... --------
yp = 0;
for j = 1:n
    yp = yp + a(j) * xp^(j-1);
end

disp("=== INTERPOLAÇÃO VIA SISTEMA LINEAR ===");
disp("Matriz de Vandermonde aumentada [V | y] (antes de eliminar):");
disp([V, y]);
disp("Matriz triangular superior (aumentada):");
disp(Ab);
disp("Coeficientes a = [a0; a1; a2; ...]:");
disp(a);
disp("p(xp) com xp =");
disp(xp);
disp(yp);
