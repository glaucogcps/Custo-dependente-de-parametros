clear; clc; close all;

fprintf(' Iniciando Teste de Análise H2 para MJLS Incerto - Exemplo 1 \n');

%% 1. Definição do Sistema MJLS Incerto (Exemplo 1 do Artigo)
nx = 2; 
nw = 2; 
nz = 3; 
nu = 1;
sigma = 2; % Número de modos Markovianos
N = 4;     % Número de vértices do politopo de transição (alpha e beta)

fprintf('Parâmetros: nx=%d, nw=%d, nz=%d, nu=%d, Modos=%d, Vértices N=%d.\n', nx, nw, nz, nu, sigma, N);

A = cell(1, sigma); B = cell(1, sigma); E = cell(1, sigma);
C = cell(1, sigma); D = cell(1, sigma);

% Matrizes do Sistema (Modos 1 e 2)
A{1} = [2 2; 
        3 1];
A{2} = [1 0; 
       0.5 1];
B{1} = [2; 1];
B{2} = [0; 0];
C{1} = [1 -1; 
        1 1; 
        0 0];
C{2} = [1 0; 
        0 1; 
        0 0];
E{1} = [0.5 0; 
        0 0.4];
E{2} = [1 0; 
        0 0.8];
D{1} = [0; 0; 1];
D{2} = [0; 0; 1];

Ez = cell(1, sigma);
Ez{1} = zeros(nz, nw); % Matriz 3x2 de zeros
Ez{2} = zeros(nz, nw); % Matriz 3x2 de zeros

% Vértices do Politopo da Matriz de Transição (Gamma)
% P(alpha, beta) = [0.9-alpha, 0.1+alpha; 0.8-beta, 0.2+beta]
% alpha in [0, 0.5], beta in [0, 0.5]
Gamma = cell(1, N);
Gamma{1} = [0.9, 0.1; 
            0.8, 0.2]; % alpha=0,   beta=0
Gamma{2} = [0.4, 0.6; 
            0.8, 0.2]; % alpha=0.5, beta=0
Gamma{3} = [0.9, 0.1; 
            0.3, 0.7]; % alpha=0,   beta=0.5
Gamma{4} = [0.4, 0.6; 
            0.3, 0.7]; % alpha=0.5, beta=0.5

% Distribuição Inicial
mu = [0 1];

disp('Sistema MJLS gerado com sucesso (Exemplo 3).');
fprintf('\n');

%% 2. Configuração das Opções do Solver
opt = struct(); 
opt.solver = 'sedumi'; % Pode usar 'mosek' se preferir
opt.op = 0;            % 0 = integral/média (rho polinomial)
opt.h2 = 0;            % Minimizar a norma H2
opt.tolerance = 1e-7;
opt.indep = 1;
vetor_cor = ['r'; 'b'; 'g'; 'c'; 'm'; 'k'];

%% 3. Teste: Variando o grau de P (deg) e degrho
if N == 2
    figure('Name', 'Análise H_2 MJLS: Custo Garantido vs Real', 'Color', 'w'); 
    hold on; grid on;
    title('Análise H_2 MJLS (Exemplo 3): Variação do grau de \rho(\alpha)'); 
    ylabel('Norma H_2'); 
    xlabel('\alpha_1');
end

H_table = [];
K_list = {}; % Inicializa lista 
disp('Iniciando Teste de Análise (Variando deg e degrho)...');

graus_rho_teste = 0:4; % Testando graus de 0 a 3 para a variável de custo rho
graus_P_teste = 1:3;   % Testando graus de 1 a 2 para a matriz de Lyapunov P

for deg = graus_P_teste
    for idx = 1:length(graus_rho_teste)
        d = graus_rho_teste(idx);
        opt.deg = deg;        % Grau de P_i(alpha)
        opt.degrho = d;       % Grau de rho(alpha)
        
        fprintf('Resolvendo para deg(P)=%d, degrho=%d... ', deg, d);
        
        out =  h2_lmi_mjls_d_incerto_controle(A, B, E, C, D, Gamma, mu, opt);
        
        if out.feas == 1
            fprintf('Viável! (Tempo: %.2f s, H2 Pior Caso: %.4f)\n', out.cpusec, out.h2);
            fprintf('Matriz de Ganho K_wc (deg=%d, d=%d):\n', deg, d);
            disp(out.K); 
            % Armazena cada modo em uma coluna da mesma linha 
            % O uso de (end+1, 1:sigma) garante que os ganhos fiquem lado a lado
            K_list(end+1, 1:sigma) = out.K(:)'; 
            
            % Plota os resultados APENAS se N == 2 e fixando um grau de P para clareza no gráfico
            if N == 2 && deg == 3
                cor_atual = vetor_cor(idx);
                % Extração dos dados padronizados
                alpha_1 = out.alpha(:, 1);
                custo_garantido = out.gcosts;
                custo_real = out.realCosts;
                % Plota o Custo REAL (linha contínua)
                plot(alpha_1, custo_real, 'Color', cor_atual, 'LineStyle', '-', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d_{\\rho} = %d', d));
                % Plota o Custo GARANTIDO (linha tracejada)
                plot(alpha_1, custo_garantido, 'Color', cor_atual, 'LineStyle', '--', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d_{\\rho} = %d', d));
            end
            
            % Calcula métricas de erro (Gap entre garantido e real)
            custo_garantido = out.gcosts;
            custo_real = out.realCosts;
            
            % Remove possíveis NaNs (caso o sistema real seja instável em algum ponto)
            idx_validos = ~isnan(custo_real);
            erro_norma = norm(custo_garantido(idx_validos) - custo_real(idx_validos));
            max_gap = max(abs(custo_garantido(idx_validos) - custo_real(idx_validos)));
            min_gap = min(abs(custo_garantido(idx_validos) - custo_real(idx_validos)));
            norma_melhor_caso = min(custo_real(idx_validos));
            norma_pior_caso = max(custo_real(idx_validos));
            
            % Armazena os dados na tabela
            H_table = [H_table; deg, d, erro_norma, max_gap, min_gap, norma_melhor_caso, norma_pior_caso, out.V, out.L];
            
            % Limpa a base de dados interna do YALMIP para evitar lentidão extrema
            % nas próximas iterações do laço (deg e degrho)
            yalmip('clear'); 
        else
            fprintf('Inviável.\n');
        end
    end
end

if N == 2
    legend('Location', 'best'); 
    xlim([0 1]);
end

%% 4. Exibição dos Resultados
disp(' ');
disp(' Tabela de Resultados (deg e degrho variando):');
if isempty(H_table)
    warning('Nenhuma solução viável foi encontrada para os graus testados.'); 
else
    T_H = array2table(H_table, 'VariableNames', {'Grau_P', 'Grau_rho', 'Erro_Norma', 'Max_Gap', 'Min_Gap', 'Norma_Melhor_Caso', 'Norma_Pior_Caso', 'Variaveis', 'Linhas_LMI'});
     % Distribui os ganhos na tabela automaticamente 
    % Como K_list já é M x sigma, basta atribuir cada coluna
    for m = 1:sigma
        col_name = sprintf('K_wc_Mode%d', m);
        T_H.(col_name) = K_list(:, m);
    end
    
    disp(T_H);
    
    writetable(T_H, 'teste_controle_H2_MJLS_Ex1_Tabela.csv'); 
%     writetable(T_H, 'teste_controle_H2_MJLS_Ex1_N2_Tabela.csv'); % vértices 1 e 4
    fprintf('Tabela salva como teste_controle_H2_MJLS_Ex1_Tabela.csv\n');
end

if N == 2 && ~isempty(H_table)
    print('teste_controle_H2_MJLS_Ex1_grafico', '-depsc'); 
    savefig('teste_controle_H2_MJLS_Ex1_grafico.fig');
%     print('teste_controle_H2_MJLS_Ex1_N2_grafico', '-depsc'); 
%     savefig('teste_controle_H2_MJLS_Ex1_N2_grafico.fig');
    fprintf('Gráfico salvo como teste_controle_H2_MJLS_Ex1_grafico.eps e teste_controle_H2_MJLS_Ex1_grafico.fig\n');
end

fprintf('Teste concluído.\n');