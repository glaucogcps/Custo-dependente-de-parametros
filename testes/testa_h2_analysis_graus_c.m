% Script de Teste para Norma H2 Contínua 
close all; clear; clc;

% 1. Carrega o arquivo de sistemas contínuos gerado anteriormente
arquivo = 'polytop_hinf_finite_grid_vs_lmi_N2_n3_step0.01_cont.mat';
if ~isfile(arquivo), error('Arquivo %s não encontrado.', arquivo); end
load(arquivo);

% 2. Seleciona um sistema específico para o teste (ex: o primeiro da lista)
idx_sys = 9; 
fprintf('Testando Sistema %d...\n', idx_sys);

% Extrai as matrizes do cell array 'matrix_configs'
Ai = matrix_configs{idx_sys,1};
Bi = matrix_configs{idx_sys,2};
Ci = matrix_configs{idx_sys,3};
Di = matrix_configs{idx_sys,4};

% Converte para double e zera D (Crítico para H2)
prec = 100000;
A = cell(1,2); B = cell(1,2); C = cell(1,2); D = cell(1,2);
for j=1:2
    A{j} = double(int32(prec*Ai{j}))/prec;
    B{j} = double(int32(prec*Bi{j}))/prec;
    C{j} = double(int32(prec*Ci{j}))/prec;
    D{j} = zeros(size(Di{j})); % D deve ser nula para H2 finita
end

% Configurações Gerais
opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; % 1 = Relaxada (Finsler), 0 = Clássica
opt.op = 0;       % 0 = Otimiza Integral (Polinomial), 1 = Otimiza Pico

%  Teste 1: Variando o grau de W (degGamma) com P fixo (deg=1) 
figure('Color', 'w'); hold on; grid on;
title(sprintf('H2 Contínua: Variação de W (P grau 1) - Sistema %d', idx_sys));
xlabel('\alpha'); ylabel('Norma H2');

H = [];
colors = jet(11); % Cores para o loop

fprintf('\n Iniciando Teste 1 (Var W, P=1) \n');
fprintf('%-5s | %-10s %-10s %-10s\n', 'Grau', 'Erro', 'Vars', 'Linhas');

for d = 0:10
    opt.deg = 1;        
    opt.degGamma = d;   
    
    try
        % Chamada da função CONTÍNUA
        out = h2_lmi_c_incerto_param(A, B, C, D, opt);
    catch
        out.feas = 0;
    end
    
    if out.feas
        % Plota Real (apenas na primeira iteração)
        if d == 0
            plot(0:0.05:1, out.realCosts, 'k--', 'LineWidth', 2, 'DisplayName', 'Real');
        end
        
        % Plota Garantido
        plot(0:0.05:1, out.gcosts, 'Color', colors(d+1,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('W grau %d', d));
         
        % Calcula erro (Norma da diferença entre vetores)
        erro = norm(out.gcosts - out.realCosts);
        H = [H; d, erro, out.V, out.L];
        
        fprintf('%-5d | %-10.4f %-10d %-10d\n', d, erro, out.V, out.L);
    else
        fprintf('%-5d | %-10s %-10s %-10s\n', d, 'Inviável', '-', '-');
    end
end
legend('show', 'Location', 'best');


%  Teste 2: Variando o grau de P (deg) com W fixo (degGamma=0) 
figure('Color', 'w'); hold on; grid on;
title(sprintf('H2 Contínua: Variação de P (W grau 0) - Sistema %d', idx_sys));
xlabel('\alpha'); ylabel('Norma H2');

Hp = [];
colors_p = lines(5);

fprintf('\n--- Iniciando Teste 2 (Var P, W=0) \n');
fprintf('%-5s | %-10s %-10s %-10s\n', 'Grau', 'Erro', 'Vars', 'Linhas');

for d = 1:5
    opt.deg = d;
    opt.degGamma = 0;
    
    try
        out = h2_lmi_c_incerto_param(A, B, C, D, opt);
    catch
        out.feas = 0;
    end
    
    if out.feas
        if d == 1
            plot(0:0.05:1, out.realCosts, 'k--', 'LineWidth', 2, 'DisplayName', 'Real');
        end
        
        plot(0:0.05:1, out.gcosts, 'Color', colors_p(d,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('P grau %d', d));
         
        erro = norm(out.gcosts - out.realCosts);
        Hp = [Hp; d, erro, out.V, out.L];
        
        fprintf('%-5d | %-10.4f %-10d %-10d\n', d, erro, out.V, out.L);
    else
        fprintf('%-5d | %-10s %-10s %-10s\n', d, 'Inviável', '-', '-');
    end
end
legend('show', 'Location', 'best');


%  Teste 3: Variando W (degGamma pares) com P fixo (deg=2) 
figure('Color', 'w'); hold on; grid on;
title(sprintf('H2 Contínua: Variação de W (P grau 2) - Sistema %d', idx_sys));
xlabel('\alpha'); ylabel('Norma H2');

Hnovo = [];
colors_n = parula(8);
idx_c = 1;

fprintf('\n--- Iniciando Teste 3 (Var W pares, P=2) \n');
fprintf('%-5s | %-10s %-10s %-10s\n', 'Grau', 'Erro', 'Vars', 'Linhas');

opt.deg = 2; % P fixo em grau 2

for d = 0:2:14
    opt.degGamma = d;
    
    try
        out = h2_lmi_c_incerto_param(A, B, C, D, opt);
    catch
        out.feas = 0;
    end
    
    if out.feas
        if d == 0
            plot(0:0.05:1, out.realCosts, 'k--', 'LineWidth', 2, 'DisplayName', 'Real');
        end
        
        plot(0:0.05:1, out.gcosts, 'Color', colors_n(idx_c,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('W grau %d', d));
         
        erro = norm(out.gcosts - out.realCosts);
        Hnovo = [Hnovo; d, erro, out.V, out.L];
        
        fprintf('%-5d | %-10.4f %-10d %-10d\n', d, erro, out.V, out.L);
    else
        fprintf('%-5d | %-10s %-10s %-10s\n', d, 'Inviável', '-', '-');
    end
    idx_c = idx_c + 1;
end
legend('show', 'Location', 'best');

%  Exibição Final das Tabelas 
fprintf('\n=== RESUMO FINAL ===\n');
disp('Tabela 1 (Var W, P=1):'); disp(array2table(H, 'VariableNames', {'Grau_W', 'Erro', 'Vars', 'Linhas'}));
disp('Tabela 2 (Var P, W=0):'); disp(array2table(Hp, 'VariableNames', {'Grau_P', 'Erro', 'Vars', 'Linhas'}));
disp('Tabela 3 (Var W pares, P=2):'); disp(array2table(Hnovo, 'VariableNames', {'Grau_W', 'Erro', 'Vars', 'Linhas'}));