clear; clc; close all;

fprintf(' Iniciando Teste de Análise H_\infty para MJLS Incerto - Exemplo 3 \n');

%% 1. Definição do Sistema MJLS Incerto (Exemplo 3 do Artigo)
nx = 1; 
nw = 1; 
nz = 1; 
sigma = 4; % Número de modos Markovianos
N = 2;     % Número de vértices do politopo de transição

fprintf('Parâmetros: nx=%d, nw=%d, nz=%d, Modos=%d, Vértices N=%d.\n', nx, nw, nz, sigma, N);

A = cell(1, sigma); E = cell(1, sigma);
Cz = cell(1, sigma); Ez = cell(1, sigma);

% Matrizes do Sistema (Modos 1 a 4 - Sistema Escalar)
A{1} = 0.999;
A{2} = 0.986;
A{3} = 1.082;
A{4} = 0.931;

for i = 1:sigma
    E{i}  = 1;
    Cz{i} = 1;
    Ez{i} = 2;
    % Adicionando matrizes de controle fictícias para o teste
    B{i} = 2*randn(); 
    D{i} = randn();
end

% Vértices do Politopo da Matriz de Transição (Gamma)
Gamma = cell(1, N);
Gamma{1} = [0.16 0.74 0.01 0.09; 
            0.59 0.02 0.11 0.28; 
            0.02 0.15 0.55 0.28;
            0.65 0.13 0.21 0.01];
            
Gamma{2} = [0.13 0.48 0.27 0.12; 
            0.48 0.35 0.08 0.09; 
            0.43 0.01 0.26 0.30;
            0.51 0.36 0.00 0.13];

% Distribuição Inicial
mu = [1 0 0 0];

disp('Sistema MJLS gerado com sucesso (Exemplo 3).');
fprintf('\n');

%% 2. Configuração das Opções do Solver
opt = struct(); 
opt.solver = 'mosek'; 
opt.op = 0;            % 0 = integral/média (rho polinomial)
opt.hinf = 0;          % Minimizar a norma Hinf
opt.tolerance = 1e-7;
vetor_cor = ['r'; 'b'; 'g'; 'c'; 'm'; 'k'];

%% 3. Teste: Variando o grau de P (deg) e degrho
if N == 2
    figure('Name', 'Análise H_\infty MJLS: Custo Garantido vs Real', 'Color', 'w'); 
    hold on; grid on;
    title('Análise H_\infty MJLS (Exemplo 3): Variação do grau de \rho(\alpha)'); 
    ylabel('Norma H_\infty'); 
    xlabel('\alpha_1');
end

H_table = [];
K_list = {}; % Inicializa lista 
disp('Iniciando Teste de Análise (Variando deg e degrho)...');

graus_rho_teste = 0:5; % Testando graus de 0 a 3 para a variável de custo rho
graus_P_teste = 1:5;   % Testando graus de 1 a 2 para a matriz de Lyapunov P

for deg = graus_P_teste
    for idx = 1:length(graus_rho_teste)
        d = graus_rho_teste(idx);
        opt.deg = deg;        % Grau de P_i(alpha)
        opt.degrho = d;       % Grau de rho(alpha)
        
        fprintf('Resolvendo para deg(P)=%d, degrho=%d... ', deg, d);
        
        % Chamada da função de análise (certifique-se de que o nome do arquivo bate com a sua função)
        out =  hinf_lmi_mjls_d_incerto_controle(A, B, E, Cz, D, Ez, Gamma, mu, opt);
        % out =  h2_lmi_mjls_d_incerto_analise(A, E, Cz, Ez, Gamma, mu, opt);
        
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
    
%     writetable(T_H, 'teste_controle_Hinf_MJLS_Ex3_Tabela.csv'); 
    writetable(T_H, 'teste_controle_Hinf_MJLS_Ex3_N2_Tabela.csv'); 
    fprintf('Tabela salva como teste_controle_Hinf_MJLS_Ex3_Tabela.csv\n');
end

if N == 2 && ~isempty(H_table)
%     print('teste_controle_Hinf_MJLS_Ex3_grafico', '-depsc'); 
%     savefig('teste_controle_Hinf_MJLS_Ex3_grafico.fig');
    print('teste_controle_Hinf_MJLS_Ex3_N2_grafico', '-depsc'); 
    savefig('teste_controle_Hinf_MJLS_Ex3_N2_grafico.fig');
    fprintf('Gráfico salvo como teste_controle_Hinf_MJLS_Ex3_grafico.eps e teste_controle_Hinf_MJLS_Ex3_grafico.fig\n');
end

fprintf('Teste concluído.\n');