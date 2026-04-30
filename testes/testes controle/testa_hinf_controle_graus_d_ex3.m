clear; clc; close all;

fprintf(' Iniciando Teste de Síntese de Controle H-infinito - Exemplo 3 \n');

%% 1. Definição do Sistema Politópico (Exemplo 3 do Artigo - 4 Vértices)
nx = 2; 
nu = 1; 
nw = 1; 
nz = 1; 
% N  = 4;
N = 2;

fprintf('Parâmetros: nx=%d, nu=%d, nw=%d, nz=%d, Vértices N=%d.\n', nx, nu, nw, nz, N);

A = cell(1, N); B = cell(1, N); E = cell(1, N);
Cz = cell(1, N); Dz = cell(1, N); Ez = cell(1, N);

% Vértice 1
A{1}  = [ 0.2970, -0.0039; 
        -0.8868, -0.8868]; 
B{1}  = [ 1; 
          0];
E{1}  = [ 0.2; 
        0.2]; 
Cz{1} = [ 0.2,  0.2]; 
Dz{1} = 1; 
Ez{1} = 0;

% Vértice 2
A{2}  = [-0.8909,  2.8938; 
        -0.2970,  0.2256]; 
B{2}  = [-3;
         -1]; 
E{2}  = [ 0.2;  
            0]; 
Cz{2} = [ 0.2,  0]; 
Dz{2} = -1; 
Ez{2} = 0;

% % Vértice 3
% A{3}  = [-0.2970,  0.8907; 
%         -0.2970,  0.8907];
% B{3}  = [-1; 
%          -1]; 
% E{3}  = [ 0; 
%           0]; 
% Cz{3} = [ 0,  0]; 
% Dz{3} = 1; 
% Ez{3} = 0;
% 
% % Vértice 4
% A{4}  = [-2.8903,  5.7876; 
%         0.9704, -6.0132]; 
% B{4}  = [-6; 
%           7]; 
% E{4}  = [ 0; 
%           0]; 
% Cz{4} = [ 0,  0]; 
% Dz{4} = -1; 
% Ez{4} = 0;

disp('Sistema gerado com sucesso (Exemplo 3 completo - 4 vértices).');
fprintf('\n');

%% 2. Configuração das Opções do Solver
opt = struct(); 
opt.solver = 'mosek'; 
opt.verbose = 0; 
opt.varFolga = 1; 
opt.op = 0;      
vetor_cor = ['r--'; 'b--'; 'g--'; 'c--'; 'm--'; 'k--'];

%% 3. Teste: Variando o grau de P (deg) de 1 a 5 e degGamma de 0 a 5
if N == 2
    figure('Name', 'Síntese H_\infty: Custo Garantido vs Real', 'Color', 'w'); 
    hold on; grid on;
    ylabel('Hinf'); 
    xlabel('\alpha');
end

H_table = [];
disp('Iniciando Teste de Síntese (Variando deg e degGamma)...');
graus_gamma_teste = 0:5;

for deg = 1:5
    for idx = 1:length(graus_gamma_teste)
        d = graus_gamma_teste(idx);
        opt.deg = deg;        % Grau de P(alpha), X(alpha), Z(alpha)
        opt.degGamma = d;     % Grau de rho(alpha)
        
        fprintf('Resolvendo para deg=%d, degGamma=%d... ', deg, d);
        out = hinf_lmi_d_incerto_controle(A, B, E, Cz, Dz, Ez, opt);
        
        if out.feas == 1
            fprintf('Viável! (Tempo: %.2f s)\n', out.cpusec_s);
            
            % Plota os resultados APENAS se N == 2
            if N == 2 && deg == 3
                % Extrai a cor atual
                cor_atual = vetor_cor(idx, 1);
                % Plota o Custo REAL usando apenas a primeira coluna de alpha (alpha_1)
                plot(out.alpha(:, 1), out.realCosts, 'Color', cor_atual, 'LineStyle', '-', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d=%d', d));
                % Plota o Custo GARANTIDO usando apenas a primeira coluna de alpha (alpha_1)
                plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('d=%d', d));
            end
            
            erro_norma = norm(out.gcosts - out.realCosts);
            max_gap = max(abs(out.gcosts - out.realCosts));
            min_gap = min(abs(out.gcosts - out.realCosts));
            norma_melhor_caso = min(out.realCosts);
            norma_pior_caso = max(out.realCosts);
            
            % Adicionando a variável 'deg' (Grau de P) como a primeira coluna da tabela
            H_table = [H_table; deg, d, erro_norma, max_gap, min_gap, norma_melhor_caso, norma_pior_caso, out.V, out.L];
        else
            fprintf('Inviável.\n');
        end
    end
end

if N == 2
    legend('Location', 'southeast'); 
    xlim([0 1]);
end

%% 4. Exibição dos Resultados
disp(' Tabela de Resultados (deg e degGamma variando):');
if isempty(H_table)
    warning('Nenhuma solução viável foi encontrada para os graus testados.'); 
    disp('O sistema não é estabilizável para os parâmetros fornecidos.');
else
    % Atualizando os nomes das colunas para incluir as novas métricas
    T_H = array2table(H_table, 'VariableNames', {'Grau_P', 'Grau_mu', 'Erro_Norma', 'Max_Gap', 'Min_Gap', 'Norma_Melhor_Caso', 'Norma_Pior_Caso', 'Variaveis', 'Linhas_LMI'});
    disp(T_H);
    
%     writetable(T_H, 'teste_controle_Hinf_Ex3_Tabela.csv'); 
    writetable(T_H, 'teste_controle_Hinf_Ex3_N2_Tabela.csv'); 
    fprintf('Tabela salva como teste_controle_Hinf_Ex3_Tabela.csv\n');
end

if N == 2 && ~isempty(H_table)
    print('teste_controle_Hinf_Ex3_N2_grafico', '-depsc'); 
    savefig('teste_controle_Hinf_Ex3_N2_grafico.fig');
    fprintf('Gráfico salvo como teste_controle_Hinf_Ex3_grafico.eps e teste_controle_Hinf_Ex3_grafico.fig\n');
elseif N ~= 2
    fprintf('Gráfico não gerado pois N = %d (plotagem suportada apenas para N=2).\n', N);
end

fprintf('Teste concluído.\n');