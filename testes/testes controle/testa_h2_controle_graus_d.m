clear; clc; close all;
fprintf(' Iniciando Teste de Síntese de Controle H2 Discreto Incerto \n');

%% 1. Geração do Sistema Politópico para Controle

nx = 3; % Número de estados
nu = 2; % Número de entradas de controle (u)
nw = 2; % Número de entradas de distúrbio (w)
nz = 2; % Número de saídas de desempenho (z)
N  = 2; % Número de vértices (N)

fprintf('Parâmetros: nx=%d, nu=%d, nw=%d, nz=%d, Vértices N=%d.\n', nx, nu, nw, nz, N);

A = cell(1, N); B = cell(1, N); E = cell(1, N);
Cz = cell(1, N); Dz = cell(1, N); Ez = cell(1, N);

for i = 1:N
    % Gera matrizes aleatórias
    A{i}  = randn(nx, nx) * 0.8; 
    B{i}  = randn(nx, nu);
    E{i}  = randn(nx, nw);
    Cz{i} = randn(nz, nx);
    Dz{i} = randn(nz, nu);
    Ez{i} = zeros(nz, nw); 
end
disp('Sistema gerado com sucesso. Matrizes Ez foram zeradas para garantir norma H2 finita.');

%% 2. Configuração das Opções do Solver
opt = struct();
opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; 
opt.op = 0;       
vetor_cor = ['r--'; 'b--'; 'g--'; 'c--'; 'm--'; 'k--'];

%% 3. Teste: Variando o grau de Gamma (degGamma) com P, H, Z fixos 

% Só cria a figura se o sistema tiver exatamente 2 vértices
if N == 2
    figure('Name', 'Síntese H2: Custo Garantido vs Real', 'Color', 'w');
    hold on; grid on;
    title('Síntese H2: Variação do grau (d) de \rho(\alpha) (Variáveis de decisão grau 3)');
    ylabel('Norma H_2 de Malha Fechada');
    xlabel('\alpha');
end

H_table = [];
K_list = {};
disp('Iniciando Teste de Síntese (Variando degGamma)...');

% Vamos testar graus de 0 a 5 para o multiplicador de custo
graus_gamma_teste = 0:5; 

for idx = 1:length(graus_gamma_teste)
    d = graus_gamma_teste(idx);   
    opt.deg = 3;        % Grau de P(alpha), H(alpha), Z(alpha)
    opt.degGamma = d;   % Grau de rho(alpha)
    
    fprintf('Resolvendo para deg=3, degGamma=%d... ', d);
    out = h2_lmi_d_incerto_controle_v2(A, B, E, Cz, Dz, Ez, opt);

    if out.feas == 1
        fprintf('Viável! (Tempo: %.2f s)\n', out.cpusec_s);
        fprintf('Matriz de Ganho K_wc (deg=%d, d=%d):\n', deg, d);
        disp(out.K_wc); 
         % Verifica se o retorno é célula (MJLS) ou matriz (Politópico)
            if iscell(out.K_wc)
                K_list(end+1, 1:numel(out.K_wc)) = out.K_wc(:)'; 
            else
                K_list{end+1, 1} = out.K_wc; % Empacota a matriz em uma célula
            end
        % Plota os resultados APENAS se N == 2
        if N == 2
            % Extrai a cor atual
            cor_atual = vetor_cor(idx, 1);
            % Plota o Custo REAL usando apenas a primeira coluna de alpha (alpha_1)
            plot(out.alpha(:, 1), out.realCosts, 'Color', cor_atual, 'LineStyle', '-', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('d=%d', d));
            % Plota o Custo GARANTIDO usando apenas a primeira coluna de alpha (alpha_1)
            plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('d=%d', d));
        end

        % Armazena dados para a tabela (Isso ocorre independente de N)
        erro_norma = norm(out.gcosts - out.realCosts);
        max_gap = max(abs(out.gcosts - out.realCosts));
        
        min_gap = min(abs(out.gcosts - out.realCosts));
        norma_melhor_caso = min(out.realCosts);
        norma_pior_caso = max(out.realCosts);

        H_table = [H_table; d, erro_norma, max_gap, min_gap, norma_melhor_caso, norma_pior_caso, out.V, out.L];
    else
        fprintf('Inviável.');
    end
end

% Formata a figura apenas se ela foi criada
if N == 2
    legend('Location', 'southeast');
    xlim([0 1]);
end

%% 4. Exibição dos Resultados

disp(' Tabela de Resultados (degGamma variando, deg=3):');
if isempty(H_table)
    warning('Nenhuma solução viável foi encontrada para os graus testados.');
    disp('O sistema gerado aleatoriamente provavelmente não é estabilizável.');
else
    % Cria a tabela apenas se houver dados
    T_H = array2table(H_table, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'Max_Gap', 'Min_Gap', 'Norma_Melhor_Caso', 'Norma_Pior_Caso', 'Variaveis', 'Linhas_LMI'});
    % Distribui os ganhos na tabela automaticamente 
    % Como K_list já é M x sigma, basta atribuir cada coluna
    for m = 1:sigma
        col_name = sprintf('K_wc_Mode%d', m);
        T_H.(col_name) = K_list(:, m);
    end
    
    disp(T_H);
    % Salva os resultados da tabela
    writetable(T_H, 'teste_controle_H2_Tabela.csv');
    fprintf('Tabela salva como teste_controle_H2_Tabela.csv\n');
end

% Salva o gráfico apenas se N == 2 e se houver dados viáveis
if N == 2 && ~isempty(H_table)
    print('teste_controle_H2_grafico', '-depsc');
    savefig('teste_controle_H2_grafico.fig'); 
    fprintf('Gráfico salvo como teste_controle_H2_grafico.eps e teste_controle_H2_grafico.fig\n');
elseif N ~= 2
    fprintf('Gráfico não gerado pois N = %d (plotagem suportada apenas para N=2).\n', N);
end

fprintf('Teste concluído.\n');