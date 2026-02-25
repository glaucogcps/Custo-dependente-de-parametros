clear; % Limpa todas as variáveis do workspace
clc;   % Limpa a janela de comando
close all; % Fecha todas as figuras abertas

fprintf(' Iniciando Geração de Politopos para Análise de Norma H2 Contínua \n');
fprintf('Buscando configurações de sistema politópico com norma H2 finita e pior caso não em vértice.\n\n');

%  1. Definição dos Parâmetros do Sistema 
n = 3; % Número de estados (ordem do sistema)
m = 2; % Número de entradas
p = 2; % Número de saídas
N = 3; % Número de vértices do politopo (quantos conjuntos de (A,B,C,D) compõem o politopo)

%  2. Definição do Passo para a Partição do Simplex 
% Um passo menor gera uma malha mais fina e mais pontos, resultando em maior precisão
% mas também maior tempo de computação.
% Para este teste, usaremos um passo de 0.01 para agilizar a busca.
% Para uma análise mais rigorosa, considere usar 0.001 ou menor.
simplex_step = 0.01; 

fprintf('Parâmetros do sistema: n=%d (estados), m=%d (entradas), p=%d (saídas), N=%d (vértices).\n', n, m, p, N);
fprintf('Usando simplex_step = %.3f para a partição do simplex.\n', simplex_step);

%  3. Configuração da Geração e Armazenamento 
matrix_configs = {}; % Inicializa um cell array vazio para armazenar as configurações válidas
target_list_size = 5; % Número desejado de configurações a serem salvas
current_list_size = 0;
iteration_count = 0; % Contador de iterações para acompanhamento do progresso

fprintf('\nIniciando busca por %d configurações de politopos que atendam aos critérios...\n', target_list_size);

% Loop para gerar e testar politopos até atingir o número desejado de configurações
while current_list_size < target_list_size
    iteration_count = iteration_count + 1;
    fprintf('Iteração %d: Tentando encontrar um politopo adequado... (Tamanho atual da lista: %d/%d)\n', iteration_count, current_list_size, target_list_size);

    % Geração de Matrizes Aleatórias para os Vértices do Politopo
    % As dimensões das matrizes são ajustadas para 'n', 'm' e 'p' definidos.
    A_cell = cell(1, N);
    B_cell = cell(1, N);
    C_cell = cell(1, N);
    D_cell = cell(1, N); 

    for i = 1:N
        A_cell{i} = randn(n, n); % Matriz A de dimensão n x n
        B_cell{i} = randn(n, m); % Matriz B de dimensão n x m
        C_cell{i} = randn(p, n); % Matriz C de dimensão p x n
        D_cell{i} = zeros(p, m); % Matriz D DEVE ser zero para a norma H2 contínua ser finita
    end

    % Chamada da Função Principal para o Caso Contínuo
    % O tipo de análise é 'continuous' e 'D_cell' é o 4º argumento.
    type_analysis_cont = 'continuous';
    raizes = lugar_raizes_matriz(A_cell, type_analysis_cont, 'simplex_step', simplex_step, 'plot_results', false);
    if raizes.max_real_part > 0
        min_val = 0.05;
        max_val = 0.5;
        r = min_val + (max_val - min_val) * rand();
        for i = 1:N
            A_cell{i} = A_cell{i} - eye(n)*(raizes.max_real_part + r);
        end
    end
    
    output_cont = norma_h2_sistema_incerto(A_cell, B_cell, C_cell, D_cell, type_analysis_cont, 'simplex_step', simplex_step, 'plot_results', false);

    %  Verificação das Condições para Salvar a Configuração 
    
    % 1. A Norma H2 de pior caso deve ser finita.
    is_finite_norm = isfinite(output_cont.max_h2_norm);

    % 2. O ponto de pior caso (alpha_star) não deve ser um vértice.
    % Um ponto é considerado um vértice se uma das componentes de alpha_star é 1 (e as outras 0).
    tol = 1e-7;
    is_not_at_vertex = ~any(abs(output_cont.alpha_star - 1) < tol);

    if is_not_at_vertex
        current_list_size = current_list_size + 1;
        % Salva as matrizes A, B, C, D e os resultados da análise em uma estrutura
        matrix_configs{current_list_size} = struct(...
            'A_cell', A_cell, ...
            'B_cell', B_cell, ...
            'C_cell', C_cell, ...
            'D_cell', D_cell, ...
            'max_h2_norm', output_cont.max_h2_norm, ...
            'alpha_star', output_cont.alpha_star);
        fprintf('  -> Politopo salvo! Tamanho da lista: %d/%d\n', current_list_size, target_list_size);
    else
        fprintf('  -> Condições não atendidas (Norma Finita: %d, Não em Vértice: %d). Tentando novamente...\n', is_finite_norm, is_not_at_vertex);
    end
end

%  4. Salvar a Lista de Configurações 
% O nome do arquivo reflete o número de vértices (N), ordem do sistema (n) e o passo do simplex.
filename = sprintf('polytop_h2_finite_non_vertex_N%d_n%d_step%g.mat', N, n, simplex_step);
save(filename, 'matrix_configs', 'N', 'n', 'm', 'p', 'simplex_step');

fprintf('\n--- Geração e Filtragem de Politopos Concluída \n');
fprintf('Salvas %d configurações de politopos em "%s"\n', current_list_size, filename);
fprintf('Parâmetros do arquivo: N=%d (vértices), n=%d (ordem do sistema), simplex_step=%g.\n', N, n, simplex_step);