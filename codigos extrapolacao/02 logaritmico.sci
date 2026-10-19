// ============================================================
// 02 - EXTRAPOLAÇÃO LOGARÍTMICA (AJUSTE LOGARÍTMICO)
//      Modelo: y = a + b * ln(x)
//      Linearização: t = ln(x)  ->  y = A + B*t  (A = a, B = b)
//      A reta é ajustada por Mínimos Quadrados: as equações
//      normais 2x2 são resolvidas por Eliminação de Gauss
//      (triangular superior + retrosubstituição).
//      Depois avalia-se o modelo em xp (extrapolação).
// ============================================================

clear; clc;

// Pontos tabelados (x deve ser > 0 para aplicar ln)
x = [1; 2; 3; 4; 5; 6; 7; 8];
y = [2; 4.079; 5.296; 6.159; 6.828; 7.375; 7.838; 8.238];

// Ponto fora do intervalo dos dados (extrapolação)
xp = 10;

// Número de pontos
n = length(x);

// -------- Linearização: t = ln(x) --------
t = log(x);

// -------- Somas para as equações normais de y = A + B*t --------
St  = sum(t);
Sy  = sum(y);
Stt = sum(t .* t);
Sty = sum(t .* y);

// -------- Equações normais ampliadas [N | s] --------
// [ n   St ] [A]   [Sy ]
// [ St  Stt] [B] = [Sty]
Ab = [n St Sy; St Stt Sty];

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

// -------- Modelo logarítmico: a = A, b = B (sem transformação) --------
a = A;
b = B;

// -------- Extrapolação: avalia o modelo em xp --------
yp = a + b * log(xp);

// Valores ajustados nos pontos tabelados (para conferir o ajuste)
y_ajustado = a + b * log(x);

disp("=== EXTRAPOLAÇÃO LOGARÍTMICA: y = a + b * ln(x) ===");
disp("Equações normais ampliadas [N | s]:");
disp([n St Sy; St Stt Sty]);
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
