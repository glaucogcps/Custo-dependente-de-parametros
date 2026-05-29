clear; clc;

%% 1. Definição do Sistema MJLS Incerto (Exemplo 3 do Artigo)
nx = 4; 
nu = 1;
nw = 2; 
nz = 4; 
sigma = 2; % Número de modos Markovianos
N = 2;     % Número de vértices do politopo de transição
fprintf('Parâmetros: nx=%d, nw=%d, nz=%d, Modos=%d, Vértices N=%d.\n', nx, nw, nz, sigma, N);
A = cell(1, sigma); B = cell(1, sigma); E = cell(1, sigma); 
Cz = cell(1, sigma); Dz = cell(1, sigma); Ez = cell(1, sigma);
% Matrizes do Sistema 

T = 0.5;
for hh=1:2
    A{hh} =eye(4) + [-0.0366, 0.0271, 0.0188, -0.4555;
                     0.0482, -1.0100, 0.0024, -4.0208;
                     0.1002, 0.3681, -0.7070, 1.4200;
                     0, 0, 1, 0]*T;
    Cz{hh} = [0, 0, 0, 10];
    Dz{hh} = 0;
    Ez{hh} = zeros(nz, nw);
end
% Vértices do Politopo da Matriz de Transição (Gamma)
%  theta1*theta2 (1-theta1)theta2  theta1(1-theta2)  (1-theta1)*(1-theta2) 
Gamma = cell(1, N);

pmin = 0.1;
pmax = 0.4;
Gamma{1} = ones(2,1)*[1-pmax, pmax];
Gamma{2} = ones(2,1)*[1-pmin, pmin];

            
B{1}=[0.4422, 3.5446, -5.5200, 0]'*T;
E{1} = B{1} * 0.2;
B{2}=[0, 0, 0, 0.000]'*T;
E{2} = B{2} * 0.2;

% Distribuição Inicial
mu = [0.5 0.5];
disp('Sistema MJLS gerado com sucesso (Exemplo 3).');
fprintf('\n');
%% 2. Configuração das Opções do Solver
opt = struct(); 
opt.solver = 'mosek'; 
opt.op = 0;            % 0 = integral/média (rho polinomial)
opt.h2 = 0;          % Minimizar a norma h2
opt.tolerance = 1e-7;
opt.indep = 0;
vetor_cor = ['r'; 'b'; 'g'; 'c'; 'm'; 'k'];

%% 3. Teste: Variando o grau de P (deg) e degrho
if N == 2
    figure('Name', 'Análise H_\infty MJLS: Custo Garantido vs Real', 'Color', 'w'); 
    hold on; grid on;
    ylabel('H_2'); 
    xlabel('prob falha');
end

H_table = [];
Graph_Points = []; % Inicializa acumulador para exportação dos pontos do gráfico
K_list = {}; % Inicializa lista 
disp('Iniciando Teste de Análise (Variando deg e degrho)...');

graus_rho_teste = 0:5; % Testando graus de 0 a 3 para a variável de custo rho
graus_P_teste = 3:5;   % Testando graus de 1 a 2 para a matriz de Lyapunov P

for deg = graus_P_teste
    for idx = 1:length(graus_rho_teste)
        d = graus_rho_teste(idx);
        opt.deg = deg;        % Grau de P_i(alpha)
        opt.degrho = d;       % Grau de rho(alpha)
        opt.degW = d;
        
        fprintf('Resolvendo para deg(P)=%d, degrho=%d... ', deg, d);
        
        % Chamada da função de análise (certifique-se de que o nome do arquivo bate com a sua função)
        out =  h2_lmi_mjls_d_incerto_controle(A, B, E, Cz, Dz, Gamma, mu, opt);
        
        if out.feas == 1
            fprintf('Viável! (Tempo: %.2f s, h2 Pior Caso: %.4f)\n', out.cpusec, out.h2);
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
                plot(pmax.*alpha_1 + pmin.*(1-alpha_1), custo_real, 'Color', cor_atual, 'LineStyle', '-', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d = %d', d));
                % Plota o Custo GARANTIDO (linha tracejada)
                plot(pmax.*alpha_1 + pmin.*(1-alpha_1), custo_garantido, 'Color', cor_atual, 'LineStyle', '--', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d = %d', d));
                % Acumula os pontos para exportação CSV
                Graph_Points = [Graph_Points; alpha_1, repmat(d, size(alpha_1)), custo_garantido, custo_real];
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
    xlim([pmin pmax]);
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
    
    writetable(T_H, 'teste_controle_h2_MJLS_NCS_N2_Tabela.csv');
%     writetable(T_H, 'teste_controle_h2_MJLS_NCS_N2_indep_Tabela.csv');
    fprintf('Tabela salva como teste_controle_h2_MJLS_NCS_Tabela.csv\n');
end

if N == 2 && ~isempty(Graph_Points)
    T_Points = array2table(Graph_Points, 'VariableNames', {'Alpha', 'Grau_rho', 'Custo_Garantido', 'Custo_Real'});
    writetable(T_Points, 'teste_controle_h2_MJLS_NCS_N2_PontosGrafico.csv');
%     writetable(T_Points, 'teste_controle_h2_MJLS_NCS_N2_indep_PontosGrafico.csv');
    fprintf('Pontos do gráfico salvos como teste_controle_H2_MJLS_NCS_N2_PontosGrafico.csv\n');
end

if N == 2 && ~isempty(H_table)
    print('teste_controle_h2_MJLS_NCS_N2_grafico', '-depsc'); 
    savefig('teste_controle_h2_MJLS_NCS_N2_grafico.fig');
%     print('teste_controle_h2_MJLS_NCS_N2_indep_grafico', '-depsc'); 
%     savefig('teste_controle_h2_MJLS_NCS_N2_indep_grafico.fig');
    fprintf('Gráfico salvo como teste_controle_h2_MJLS_Ex3_grafico.eps e teste_controle_h2_MJLS_Ex3_grafico.fig\n');
end

fprintf('Teste concluído.\n');