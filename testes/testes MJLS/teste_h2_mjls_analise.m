% % Exemplo 3 - 
% https://www.tandfonline.com/doi/epdf/10.1080/00207170802136178?needAccess=true
% clear; clc; close all;
% A{1} = 0.999;
% A{2} = 0.986;
% A{3} = 1.082;
% A{4} = 0.931;
% for i=1:4
% E{i} = 1;
% Cz{i} = 1;
% Ez{i} = 0;
% end
% Gamma{1} = [0.16 0.74 0.01 0.09; 
%             0.59 0.02 0.11 0.28; 
%             0.02 0.15 0.55 0.28;
%             0.65 0.13 0.21 0.01];
% Gamma{2} = [0.13 0.48 0.27 0.12; 
%             0.48 0.35 0.08 0.09; 
%             0.43 0.01 0.26 0.30;
%             0.51 0.36 0.00 0.13];
% mu = [1 0 0 0];
% output=h2_lmi_mjls_d_incerto(A,E,Cz,Ez,Gamma,mu)
% if output.feas==1
% plot(output.realCosts(:,1),output.realCosts(:,end))
% hold on
% plot(output.realCosts(:,1),output.realCosts(:,2),'r--')
% end

clear; clc; close all;

fprintf(' Iniciando Teste de Análise H2 para MJLS Incerto - Exemplo 2 \n');

%% 1. Definição do Sistema MJLS Incerto (Exemplo 2 do Artigo)
nx = 2; 
nw = 2; 
nz = 2; 
sigma = 3; % Número de modos Markovianos
% N = 2;     % Número de vértices do politopo de transição
N = 3;

fprintf('Parâmetros: nx=%d, nw=%d, nz=%d, Modos=%d, Vértices N=%d.\n', nx, nw, nz, sigma, N);

A = cell(1, sigma); E = cell(1, sigma);
Cz = cell(1, sigma); Ez = cell(1, sigma);

% Matrizes do Sistema (Modos 1 a 3)
A{1} = [0.6 0.4; 0.3 0.7];
A{2} = [0.2 0.6; 0.6 0.5];
A{3} = [0.3 0.3; 0.4 0.5];

for i = 1:sigma
    E{i}  = [2 0; 0 1];
    Cz{i} = eye(2);
    Ez{i} = zeros(2,2);
end

% Vértices do Politopo da Matriz de Transição (Gamma)
Gamma = cell(1, N);
Gamma{1} = [0.09 0.81 0.10; 
            0.27 0.32 0.41; 
            0.11 0.32 0.57];
            
Gamma{2} = [0.74 0.04 0.22; 
            0.26 0.36 0.38; 
            0.71 0.02 0.27];

Gamma{3} = [0.64 0.33 0.03; 
            0.17 0.36 0.47; 
            0.79 0.20 0.01];

% Distribuição Inicial
mu = [1 0 0];

disp('Sistema MJLS gerado com sucesso.');
fprintf('\n');

%% 2. Configuração das Opções do Solver
opt = struct(); 
opt.solver = 'mosek'; 
opt.op = 0;            % 0 = integral/média (rho polinomial)
opt.h2 = 0;            % Minimizar a norma H2
opt.tolerance = 1e-7;
vetor_cor = ['r'; 'b'; 'g'; 'c'; 'm'; 'k'];

%% 3. Teste: Variando o grau de P (deg) e degrho
if N == 2
    figure('Name', 'Análise H_2 MJLS: Custo Garantido vs Real', 'Color', 'w'); 
    hold on; grid on;
    title('Análise H_2 MJLS: Variação do grau de \rho(\alpha)'); 
    ylabel('Norma H_2'); 
    xlabel('\alpha_1');
end

H_table = [];
disp('Iniciando Teste de Análise (Variando deg e degrho)...');

graus_rho_teste = 0:3; % Testando graus de 0 a 3 para a variável de custo rho
graus_P_teste = 1:2;   % Testando graus de 1 a 2 para a matriz de Lyapunov P

for deg = graus_P_teste
    for idx = 1:length(graus_rho_teste)
        d = graus_rho_teste(idx);
        opt.deg = deg;        % Grau de P_i(alpha)
        opt.degrho = d;       % Grau de rho(alpha)
        
        fprintf('Resolvendo para deg(P)=%d, degrho=%d... ', deg, d);
        out = h2_lmi_mjls_d_incerto_analise(A, E, Cz, Ez, Gamma, mu, opt);
        
        if out.feas == 1
            fprintf('Viável! (Tempo: %.2f s, H2 Pior Caso: %.4f)\n', out.cpusec, out.h2);
            
            % Plota os resultados APENAS se N == 2 e fixando um grau de P para clareza no gráfico
            if N == 2 && deg == 1
                cor_atual = vetor_cor(idx);
                % Extração dos dados da matriz out.realCosts
                % Formato: [alpha_1, alpha_2, Custo_Garantido, Custo_Real]
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
            erro_norma = norm(custo_garantido - custo_real);
            max_gap = max(abs(custo_garantido - custo_real));
            min_gap = min(abs(custo_garantido - custo_real));
            norma_melhor_caso = min(custo_real);
            norma_pior_caso = max(custo_real);
            
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
    disp(T_H);
    
%     writetable(T_H, 'teste_analise_H2_MJLS_Ex2_N2_Tabela.csv'); 
    writetable(T_H, 'teste_analise_H2_MJLS_Ex2_Tabela.csv');
    fprintf('Tabela salva como teste_analise_H2_MJLS_Ex2_Tabela.csv\n');
end

if N == 2 && ~isempty(H_table)
%     print('teste_analise_H2_MJLS_Ex2_N2_grafico', '-depsc'); 
%     savefig('teste_analise_H2_MJLS_Ex2_N2_grafico.fig');
    print('teste_analise_H2_MJLS_Ex2_grafico', '-depsc'); 
    savefig('teste_analise_H2_MJLS_Ex2_grafico.fig');
    fprintf('Gráfico salvo como teste_analise_H2_MJLS_Ex2_grafico.eps e teste_analise_H2_MJLS_Ex2_grafico.fig\n');
end

fprintf('Teste concluído.\n');