clear; clc; close all;

% Carrega arquivo
arquivo = 'polytop_hinf_finite_grid_vs_lmi_N2_n3_step0.01_cont.mat';
if ~isfile(arquivo), error('Arquivo não encontrado'); end
load(arquivo);

% Configurações
% 5 é bom, 9 melhor
idx_sys = 10; % Escolha o sistema
prec = 100000;
graus_teste = 1:4; % Testar graus 1, 2, 3, 4

fprintf('Gerando tabela para Sistema %d (Comparação Constante vs Polinomial)...\n', idx_sys);

% Extração e Conversão
Ai = matrix_configs{idx_sys,1}; Bi = matrix_configs{idx_sys,2};
Ci = matrix_configs{idx_sys,3}; Di = matrix_configs{idx_sys,4};
A = cell(1,2); B = cell(1,2); C = cell(1,2); D = cell(1,2);
for j=1:2
    A{j} = double(int32(prec*Ai{j}))/prec; B{j} = double(int32(prec*Bi{j}))/prec;
    C{j} = double(int32(prec*Ci{j}))/prec; D{j} = double(int32(prec*Di{j}))/prec;
end

% 1. Curva Real
alphas = 0:0.01:1;
H_real = zeros(length(alphas), 1);
for k = 1:length(alphas)
    a = alphas(k);
    sys = ss(a*A{1}+(1-a)*A{2}, a*B{1}+(1-a)*B{2}, a*C{1}+(1-a)*C{2}, a*D{1}+(1-a)*D{2}, 0);
    H_real(k) = norm(sys, inf);
end

% Prepara Tabela e Figura
tab = [];
figure('Color', 'w'); hold on; grid on;
plot(alphas, H_real, 'k--', 'LineWidth', 2, 'DisplayName', 'Real');
colors = lines(length(graus_teste));

for i = 1:length(graus_teste)
    g = graus_teste(i);
    fprintf('Calculando Grau g = %d... ', g);
    
    %  CASO 1: Custo Constante (degGamma = 0) 
    opt1 = struct('deg', g, 'degGamma', 0, 'varFolga', 1, 'verbose', 0, 'op', 0);
    try
        out = hinf_lmi_c_incerto_param(A, B, C, D, opt1);
    catch
        out.feas = 0; out.V = 0; out.L = 0;
    end
    
    if out.feas
        curve_const = sqrt(out.wc) * ones(size(alphas));
        diff_const = norm(curve_const - H_real);
    else
        diff_const = NaN;
    end
    
    %  CASO 2: Custo Polinomial de Grau g (degGamma = g) 
    opt2 = struct('deg', g, 'degGamma', g, 'varFolga', 1, 'verbose', 0, 'op', 0);
    try
        out_poly = hinf_lmi_c_incerto_param(A, B, C, D, opt2);
    catch
        out_poly.feas = 0; out_poly.V = 0; out_poly.L = 0;
    end
    
    curve_poly = zeros(length(alphas), 1);
    if out_poly.feas
        for k = 1:length(alphas)
            val_rho = double(evalpar(out_poly.rho_alpha, {[alphas(k), 1-alphas(k)]}));
            curve_poly(k) = sqrt(max(0, val_rho));
        end
        diff_poly = norm(curve_poly - H_real);
        
        % Plota curva polinomial
        plot(alphas, curve_poly, 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Polinomial g=%d', g));
    else
        diff_poly = NaN;
    end
    
    % Tabela: [Grau, Erro_Const, Vars_Const, Erro_Poly, Vars_Poly]
    tab = [tab; g, diff_const, out.V, diff_poly, out_poly.V];
    fprintf('Ok.\n');
end

xlabel('\alpha'); ylabel('H_\infty');
title(sprintf('Sistema %d: Custo Constante vs Polinomial', idx_sys));
legend('show', 'Location', 'best');

%  Exibe Tabela 
fprintf('\n%-5s | %-12s %-8s | %-12s %-8s\n', 'Grau', 'Erro(Const)', 'Vars', 'Erro(Poly)', 'Vars');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:size(tab, 1)
    fprintf('%-5d | %-12.4f %-8d | %-12.4f %-8d\n', ...
        tab(i,1), tab(i,2), tab(i,3), tab(i,4), tab(i,5));
end

% Gráfico de Barras Comparativo
figure('Name', 'Comparação de Erros', 'Color', 'w');
bar(tab(:,1), [tab(:,2), tab(:,4)]);
legend('Custo Constante (\gamma)', 'Custo Polinomial (\gamma(\alpha))');
xlabel('Grau (g)'); ylabel('Erro L2');
title('Redução do Conservadorismo');
grid on;