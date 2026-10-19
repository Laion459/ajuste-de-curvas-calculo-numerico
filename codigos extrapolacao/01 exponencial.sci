// ============================================================
// 01 - EXTRAPOLAÇÃO EXPONENCIAL (AJUSTE EXPONENCIAL)
//      Modelo: y = a * e^(b*x)
//      Linearização: ln(y) = ln(a) + b*x  ->  z = A + B*x
//      A reta z é ajustada por Mínimos Quadrados: as equações
//      normais 2x2 são resolvidas por Eliminação de Gauss
//      (triangular superior + retrosubstituição).
//      Depois avalia-se o modelo em xp (extrapolação).
// ============================================================

clear; clc;

// Pontos tabelados (y deve ser > 0 para aplicar ln)
x = [-1; -0.7; -0.4; -0.1; 0.2; 0.5; 0.8; 1];
y = [36.547; 17.264; 8.155; 3.852; 1.82; 0.86; 0.406; 0.246];

// Ponto fora do intervalo dos dados (extrapolação)
xp = 0.6;

// Número de pontos
n = length(x);

// -------- Linearização: z = ln(y) --------
z = log(y);

// -------- Somas para as equações normais de z = A + B*x --------
Sx  = sum(x);
Sz  = sum(z);
Sxx = sum(x .* x);
Sxz = sum(x .* z);

// -------- Equações normais ampliadas [N | s] --------
// [ n   Sx ] [A]   [Sz ]
// [ Sx  Sxx] [B] = [Sxz]
Ab = [n Sx Sz; Sx Sxx Sxz];

// Tamanho do sistema (2 incógnitas: A e B)
m = size(Ab, 1);

// -------- FASE 1: Eliminação para frente --------
// Transforma N em matriz triangular superior
for k = 1:m-1
    for i = k+1:m
        fator = Ab(i, k) / Ab(k, k);
        Ab(i, k:m+1) = Ab(i, k:m+1) - fator * Ab(k, k:m+1);
    end
end

// -------- FASE 2: Substituição retroativa --------
// Última incógnita (B) já está resolvida na última linha
coef(m) = Ab(m, m+1) / Ab(m, m);

// Sobe linha por linha resolvendo A
for i = m-1:-1:1
    soma = 0;
    for j = i+1:m
        soma = soma + Ab(i, j) * coef(j);
    end
    coef(i) = (Ab(i, m+1) - soma) / Ab(i, i);
end

A = coef(1);
B = coef(2);

// -------- Volta ao modelo exponencial: a = e^A, b = B --------
a = exp(A);
b = B;

// -------- Extrapolação: avalia o modelo em xp --------
yp = a * exp(b * xp);

// Valores ajustados nos pontos tabelados (para conferir o ajuste)
y_ajustado = a * exp(b * x);

disp("=== EXTRAPOLAÇÃO EXPONENCIAL: y = a * e^(b*x) ===");
disp("Equações normais ampliadas [N | s]:");
disp([n Sx Sz; Sx Sxx Sxz]);
disp("Matriz triangular superior (aumentada):");
disp(Ab);
disp("Coeficientes da reta (retrosubstituição) [A; B] =");
disp(coef);
disp("a =");
disp(a);
disp("b =");
disp(b);
disp("y ajustado nos pontos tabelados:");
disp(y_ajustado);
disp("Extrapolação em xp =");
disp(xp);
disp("yp =");
disp(yp);
