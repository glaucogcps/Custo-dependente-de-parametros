clear; clc; close all;

fprintf(' Iniciando Teste de Análise H_infty para MJLS Incerto - Exemplo 3 \n');

%% 1. Definição do Sistema MJLS Incerto (Exemplo 3 do Artigo)
nx = 2; 
nw = 1; 
nz = 1; 
sigma = 2; % Número de modos Markovianos
N = 2;     % Número de vértices do politopo de transição

fprintf('Parâmetros: nx=%d, nw=%d, nz=%d, Modos=%d, Vértices N=%d.\n', nx, nw, nz, sigma, N);

A = cell(1, sigma); E = cell(1, sigma);
Cz = cell(1, sigma); Ez = cell(1, sigma);

beta = 2.72;
% Matrizes do Sistema 
A{1} = beta*[0.32, -0.40; 0.80, -0.80];
A{2} = beta*[0.08, -0.26; 0.80, -1.12];

B{1} = [2; 1];
B{2} = [1; -1];

E{1} = [0.08; 0.10];
E{2} = E{1};

Cz{1} = [0.5, 0];
Cz{2} = Cz{1};

Dz{1} = 0.4;
Dz{2} = 0.4;

Ez{1} = 0.6;
Ez{2} = 0.6;

% Vértices do Politopo da Matriz de Transição (Gamma)
Gamma = cell(1, N);
Gamma{1} = [1, 0;
            0.3, 0.7];
            
Gamma{2} = [0, 1;
            0.3, 0.7];

% Distribuição Inicial
mu = [1 0];

disp('Sistema MJLS gerado com sucesso (Exemplo 3).');
fprintf('\n');

%% 2. Configuração das Opções do Solver
opt = struct(); 
opt.solver = 'mosek'; 
opt.op = 0;            % 0 = integral/média (rho polinomial)
opt.hinf = 0;          % Minimizar a norma Hinf
opt.tolerance = 1e-7;
opt.indep = 1;
vetor_cor = ['r'; 'b'; 'g'; 'c'; 'm'; 'k'];

%% 3. Teste: Variando o grau de P (deg) e degrho
if N == 2
    figure('Name', 'Análise H_\infty MJLS: Custo Garantido vs Real', 'Color', 'w'); 
    hold on; grid on;
    ylabel('Hinf'); 
    xlabel('\alpha');
end

H_table = [];
K_list = {}; % Inicializa lista 
disp('Iniciando Teste de Análise (Variando deg e degrho)...');

graus_rho_teste = 0:5; % Testando graus de 0 a 3 para a variável de custo rho
graus_P_teste = 0:5;   % Testando graus de 1 a 2 para a matriz de Lyapunov P

for deg = graus_P_teste
    for idx = 1:length(graus_rho_teste)
        d = graus_rho_teste(idx);
        opt.deg = deg;        % Grau de P_i(alpha)
        opt.degrho = d;       % Grau de rho(alpha)
        
        fprintf('Resolvendo para deg(P)=%d, degrho=%d... ', deg, d);
        
        % Chamada da função de análise (certifique-se de que o nome do arquivo bate com a sua função)
        out =  hinf_lmi_mjls_d_incerto_controle(A, B, E, Cz, Dz, Ez, Gamma, mu, opt);
        
        if out.feas == 1
            fprintf('Viável! (Tempo: %.2f s, Hinf Pior Caso: %.4f)\n', out.cpusec, out.hinf);
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
                    'LineWidth', 1.5, 'DisplayName', sprintf('d = %d', d));
                % Plota o Custo GARANTIDO (linha tracejada)
                plot(alpha_1, custo_garantido, 'Color', cor_atual, 'LineStyle', '--', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d = %d', d));
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
        else
            fprintf('Inviável.\n');
        end
        
        % Limpa a base de dados interna do YALMIP para evitar lentidão extrema
        % nas próximas iterações do laço (deg e degrho)
        yalmip('clear'); 
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
    
%     writetable(T_H, 'teste_controle_Hinf_MJLS_Ex7_N2_Tabela.csv');
    writetable(T_H, 'teste_controle_Hinf_MJLS_Ex7_N2_indep_Tabela.csv');
    fprintf('Tabela salva como teste_controle_Hinf_MJLS_Ex7_Tabela.csv\n');
end

if N == 2 && ~isempty(H_table)
%     print('teste_controle_Hinf_MJLS_Ex7_N2_grafico', '-depsc'); 
%     savefig('teste_controle_Hinf_MJLS_Ex7_N2_grafico.fig');
    print('teste_controle_Hinf_MJLS_Ex7_N2_indep_grafico', '-depsc'); 
    savefig('teste_controle_Hinf_MJLS_Ex7_N2_indep_grafico.fig');
    fprintf('Gráfico salvo como teste_controle_Hinf_MJLS_Ex3_grafico.eps e teste_controle_Hinf_MJLS_Ex3_grafico.fig\n');
end

fprintf('Teste concluído.\n');