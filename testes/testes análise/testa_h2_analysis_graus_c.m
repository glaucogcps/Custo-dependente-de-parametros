% Script de Teste para Norma H2 Contínua 
clear; clc; close all;

fprintf(' Iniciando Teste de Análise H2 Contínua Incerta \n');

%% 1. Carregamento do Sistema
arquivo = 'polytop_hinf_finite_grid_vs_lmi_N2_n3_step0.01_cont.mat';
if ~isfile(arquivo), error('Arquivo %s não encontrado.', arquivo); end
load(arquivo);

% Seleciona um sistema específico para o teste
idx_sys = 9; 
fprintf('Testando Sistema %d...', idx_sys);

% Extrai as matrizes do cell array 'matrix_configs'
Ai = matrix_configs{idx_sys,1};
Bi = matrix_configs{idx_sys,2};
Ci = matrix_configs{idx_sys,3};
Di = matrix_configs{idx_sys,4};

% Define N dinamicamente com base no sistema carregado
N = length(Ai);
fprintf('Parâmetros: Vértices N=%d.', N);

% Converte para double e zera D (Crítico para H2)
prec = 100000;
A = cell(1,N); B = cell(1,N); C = cell(1,N); D = cell(1,N);
for j=1:N
    A{j} = double(int32(prec*Ai{j}))/prec;
    B{j} = double(int32(prec*Bi{j}))/prec;
    C{j} = double(int32(prec*Ci{j}))/prec;
    D{j} = zeros(size(Di{j})); % D deve ser nula para H2 finita
end
disp('Sistema carregado com sucesso. Matrizes D foram zeradas para garantir norma H2 finita.');

%% 2. Configuração das Opções do Solver
opt = struct();
opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; 
opt.op = 0;      

%% 3. Teste 1: Variando o grau de W (degGamma) com P fixo (deg=1)
fprintf('Iniciando Teste 1 (Var W, P=1)');
if N == 2
    fig1 = figure('Name', 'H2 Contínua: Teste 1', 'Color', 'w');
    hold on; grid on;
    title(sprintf('H2 Contínua: Variação de W (P grau 1) - Sistema %d', idx_sys));
    ylabel('Norma H_2'); xlabel('\alpha_1');
end

H1_table = [];
graus_teste1 = 0:10;
colors1 = jet(length(graus_teste1));

for idx = 1:length(graus_teste1)
    d = graus_teste1(idx);
    opt.deg = 1;        
    opt.degGamma = d;   
    
    fprintf('Resolvendo para deg=1, degGamma=%d... ', d);
    out = h2_lmi_c_incerto_param(A, B, C, D, opt);
    
    if out.feas == 1
        fprintf('Viável! (Tempo: %.2f s)\n', out.cpusec_s);
        
        if N == 2
            cor_atual = colors1(idx, :);
            plot(out.alpha(:, 1), out.realCosts, 'Color', cor_atual, 'LineStyle', '-', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Real (W=%d)', d));
            plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Garantido (W=%d)', d));
        end
        
        erro_norma = norm(out.gcosts - out.realCosts);
        max_gap = max(abs(out.gcosts - out.realCosts));
        H1_table = [H1_table; d, erro_norma, max_gap, out.V, out.L];
    else
        fprintf('Inviável.\n');
    end
end
if N == 2
    legend('Location', 'best'); xlim([0 1]);
end

%% 4. Teste 2: Variando o grau de P (deg) com W fixo (degGamma=0)
fprintf('Iniciando Teste 2 (Var P, W=0)');
if N == 2
    fig2 = figure('Name', 'H2 Contínua: Teste 2', 'Color', 'w');
    hold on; grid on;
    title(sprintf('H2 Contínua: Variação de P (W grau 0) - Sistema %d', idx_sys));
    ylabel('Norma H_2'); xlabel('\alpha_1');
end

H2_table = [];
graus_teste2 = 1:5;
colors2 = lines(length(graus_teste2));

for idx = 1:length(graus_teste2)
    d = graus_teste2(idx);
    opt.deg = d;        
    opt.degGamma = 0;   
    
    fprintf('Resolvendo para deg=%d, degGamma=0... ', d);
    try
        out = h2_lmi_c_incerto_param(A, B, C, D, opt);
    catch
        out.feas = 0;
    end
    
    if out.feas == 1
        fprintf('Viável! (Tempo: %.2f s)', out.cpusec_s);
        
        if N == 2
            cor_atual = colors2(idx, :);
            plot(out.alpha(:, 1), out.realCosts, 'Color', cor_atual, 'LineStyle', '-', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Real (P=%d)', d));
            plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Garantido (P=%d)', d));
        end
        
        erro_norma = norm(out.gcosts - out.realCosts);
        max_gap = max(abs(out.gcosts - out.realCosts));
        H2_table = [H2_table; d, erro_norma, max_gap, out.V, out.L];
    else
        fprintf('Inviável.');
    end
end
if N == 2
    legend('Location', 'best'); xlim([0 1]);
end

%% 5. Teste 3: Variando W (degGamma) com P fixo (deg=2)
fprintf('Iniciando Teste 3 (Var W, P=2)');
if N == 2
    fig3 = figure('Name', 'H2 Contínua: Teste 3', 'Color', 'w');
    hold on; grid on;
    title(sprintf('H2 Contínua: Variação de W (P grau 2) - Sistema %d', idx_sys));
    ylabel('Norma H_2'); xlabel('\alpha_1');
end

H3_table = [];
graus_teste3 = 0:10;
colors3 = parula(length(graus_teste3));

for idx = 1:length(graus_teste3)
    d = graus_teste3(idx);
    opt.deg = 2;        
    opt.degGamma = d;   
    
    fprintf('Resolvendo para deg=2, degGamma=%d... ', d);
    out = h2_lmi_c_incerto_param(A, B, C, D, opt);
    
    if out.feas == 1
        fprintf('Viável! (Tempo: %.2f s)', out.cpusec_s);
        
        if N == 2
            cor_atual = colors3(idx, :);
            plot(out.alpha(:, 1), out.realCosts, 'Color', cor_atual, 'LineStyle', '-', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Real (W=%d)', d));
            plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Garantido (W=%d)', d));
        end
        
        erro_norma = norm(out.gcosts - out.realCosts);
        max_gap = max(abs(out.gcosts - out.realCosts));
        H3_table = [H3_table; d, erro_norma, max_gap, out.V, out.L];
    else
        fprintf('Inviável.');
    end
end
if N == 2
    legend('Location', 'best'); xlim([0 1]);
end

%% 6. Exibição dos Resultados
fprintf(' RESUMO FINAL DOS RESULTADOS');

% Tabela 1
disp('Tabela 1 (Var W, P=1)');
if isempty(H1_table)
    disp('Nenhuma solução viável encontrada para o Teste 1.');
else
    T1 = array2table(H1_table, 'VariableNames', {'Grau_W', 'Erro_Norma', 'Max_Gap', 'Variaveis', 'Linhas_LMI'});
    disp(T1);
    writetable(T1, 'teste_H2_cont_Tabela1.csv');
end

% Tabela 2
disp('Tabela 2 (Var P, W=0)');
if isempty(H2_table)
    disp('Nenhuma solução viável encontrada para o Teste 2.');
else
    T2 = array2table(H2_table, 'VariableNames', {'Grau_P', 'Erro_Norma', 'Max_Gap', 'Variaveis', 'Linhas_LMI'});
    disp(T2);
    writetable(T2, 'teste_H2_cont_Tabela2.csv');
end

% Tabela 3
disp('Tabela 3 (Var W, P=2)');
if isempty(H3_table)
    disp('Nenhuma solução viável encontrada para o Teste 3.');
else
    T3 = array2table(H3_table, 'VariableNames', {'Grau_W', 'Erro_Norma', 'Max_Gap', 'Variaveis', 'Linhas_LMI'});
    disp(T3);
    writetable(T3, 'teste_H2_cont_Tabela3.csv');
end

% Salva os gráficos
if N == 2
    if ~isempty(H1_table), figure(fig1); print('teste_H2_cont_grafico1', '-depsc'); end
    if ~isempty(H2_table), figure(fig2); print('teste_H2_cont_grafico2', '-depsc'); end
    if ~isempty(H3_table), figure(fig3); print('teste_H2_cont_grafico3', '-depsc'); end
    fprintf('Gráficos e tabelas salvos no diretório atual.');
else
    fprintf('Gráficos não gerados pois N = %d (plotagem suportada apenas para N=2).', N);
end
fprintf('Teste concluído.');