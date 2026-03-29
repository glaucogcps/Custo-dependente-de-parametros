% Script de Teste para Análise Hinf Discreta
clear; clc; close all;

fprintf(' Iniciando Teste de Análise Hinf Discreta Incerta');

%% 1. Carregamento do Sistema
load sistema_teste_d.mat;

% Define N dinamicamente com base no sistema carregado
N = length(A);
fprintf('Parâmetros: Vértices N=%d.\n', N);
disp('Sistema carregado com sucesso.');

%% 2. Configuração das Opções do Solver
opt = struct();
opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; % Usando a formulação relaxada (Finsler)

%% 3. Teste 1: Variando o grau de Gamma (degGamma) com P fixo (deg=1) 
fprintf('Iniciando Teste 1 (Var mu, P=1)');
if N == 2
    fig1 = figure('Name', 'Hinf Discreta: Teste 1', 'Color', 'w');
    hold on; grid on;
    title('Variação do grau (d) de \mu(\alpha) (P grau 1)');
    ylabel('Norma H_\infty'); xlabel('\alpha_1');
end

H1_table = [];
graus_teste1 = 0:10;
colors1 = jet(length(graus_teste1));

for idx = 1:length(graus_teste1)
    d = graus_teste1(idx);
    opt.deg = 1;        
    opt.degGamma = d;   
    
    fprintf('Resolvendo para deg=1, degGamma=%d... ', d);
    out = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    
    if out.feas == 1
        fprintf('Viável! (Tempo: %.2f s)', out.cpusec_s);
        
        if N == 2
            cor_atual = colors1(idx, :);
            if idx == 1
                plot(out.alpha(:, 1), out.realCosts, 'k-', 'LineWidth', 2, 'DisplayName', 'Real');
            end
            plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Garantido (d=%d)', d));
        end
        
        erro_norma = norm(out.gcosts - out.realCosts);
        max_gap = max(abs(out.gcosts - out.realCosts));
        H1_table = [H1_table; d, erro_norma, max_gap, out.V, out.L];
    else
        fprintf('Inviável.');
    end
end
if N == 2
    legend('Location', 'best'); xlim([0 1]);
end

%% 4. Teste 2: Variando o grau de P (deg) comparando Gamma fixo (0, 1 e 8) 
fprintf('Iniciando Teste 2 (Var P, mu = 0, 1 e 8)');
if N == 2
    fig2 = figure('Name', 'Hinf Discreta: Teste 2', 'Color', 'w');
    hold on; grid on;
    title('Variação do grau (g) de P (Plotando \mu=0)');
    ylabel('Norma H_\infty'); xlabel('\alpha_1');
end

H2_table = [];
graus_teste2 = 1:5;
colors2 = lines(length(graus_teste2));

for idx = 1:length(graus_teste2)
    d = graus_teste2(idx);
    fprintf('Resolvendo para deg=%d...', d);
    
    try
        % Caso 1: degGamma = 0
        opt.deg = d; opt.degGamma = 0;
        out1 = hinf_lmi_d_incerto_param(A, B, C, D, opt);
        
        % Caso 2: degGamma = 1 
        opt.degGamma = 1;
        out2 = hinf_lmi_d_incerto_param(A, B, C, D, opt);
        
        % Caso 3: degGamma = 8
        opt.degGamma = 8;
        out3 = hinf_lmi_d_incerto_param(A, B, C, D, opt);
        
        if out1.feas && out2.feas && out3.feas
            fprintf('  -> Viável para todos os casos de mu');
            
            erro1 = norm(out1.gcosts - out1.realCosts); max_gap1 = max(abs(out1.gcosts - out1.realCosts));
            erro2 = norm(out2.gcosts - out2.realCosts); max_gap2 = max(abs(out2.gcosts - out2.realCosts));
            erro3 = norm(out3.gcosts - out3.realCosts); max_gap3 = max(abs(out3.gcosts - out3.realCosts));
            
            if N == 2
                cor_atual = colors2(idx, :);
                if idx == 1
                    plot(out1.alpha(:, 1), out1.realCosts, 'k-', 'LineWidth', 2, 'DisplayName', 'Real');
                end
                % Plota apenas o Caso 1 (mu=0) para não poluir o gráfico
                plot(out1.alpha(:, 1), out1.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('Garantido (g=%d, mu=0)', d));
            end
            
            H2_table = [H2_table; d, erro1, max_gap1, out1.V, out1.L, ...
                                     erro2, max_gap2, out2.V, out2.L, ...
                                     erro3, max_gap3, out3.V, out3.L];
        else
            fprintf('Inviável em pelo menos um dos casos.');
        end
    catch
        fprintf('Erro numérico ou falta de memória ao montar as LMIs.');
    end
end
if N == 2
    legend('Location', 'best'); xlim([0 1]);
end

%% 5. Gráfico de Barras Agrupadas (Teste 2)
if ~isempty(H2_table)
    fig_bar = figure('Name', 'Hinf Discreta: Comparação de Erros', 'Color', 'w');
    grau = H2_table(:,1);
    diferenca_caso1 = H2_table(:,2);  % Erro com degGamma = 0
    diferenca_caso2 = H2_table(:,6);  % Erro com degGamma = 1
    diferenca_caso3 = H2_table(:,10); % Erro com degGamma = 8
    
    bar(grau, [diferenca_caso1, diferenca_caso2, diferenca_caso3]);
    xlabel('Grau de P (g)');
    ylabel('Norma da diferença (Custo Garantido - H_\infty)');
    title('Comparação H_\infty: Diferença vs Grau de P');
    legend('\mu(\alpha) grau 0', '\mu(\alpha) grau 1', '\mu(\alpha) grau 8', 'Location', 'best');
    grid on;
end

%% 6. Teste 3: Variando Gamma (degGamma) com P fixo (deg=2) 
fprintf('Iniciando Teste 3 (Var mu, P=2)');
if N == 2
    fig3 = figure('Name', 'Hinf Discreta: Teste 3', 'Color', 'w');
    hold on; grid on;
    title('Variação do grau (d) de \mu(\alpha) (P grau 2)');
    ylabel('Norma H_\infty'); xlabel('\alpha_1');
end

H3_table = [];
graus_teste3 = 0:8;
colors3 = parula(length(graus_teste3));

for idx = 1:length(graus_teste3)
    d = graus_teste3(idx);
    opt.deg = 2;        
    opt.degGamma = d;   
    
    fprintf('Resolvendo para deg=2, degGamma=%d... ', d);
    out = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    
    if out.feas == 1
        fprintf('Viável! (Tempo: %.2f s)', out.cpusec_s);
        
        if N == 2
            cor_atual = colors3(idx, :);
            if idx == 1
                plot(out.alpha(:, 1), out.realCosts, 'k-', 'LineWidth', 2, 'DisplayName', 'Real');
            end
            plot(out.alpha(:, 1), out.gcosts, 'Color', cor_atual, 'LineStyle', '--', ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Garantido (d=%d)', d));
        end
        
        erro_norma = norm(out.gcosts - out.realCosts);
        max_gap = max(abs(out.gcosts - out.realCosts));
        H3_table = [H3_table; d, erro_norma, max_gap, out.V, out.L];
    else
        fprintf('Inviável.\n');
    end
end
if N == 2
    legend('Location', 'best'); xlim([0 1]);
end

%% 7. Exibição dos Resultados
fprintf(' RESUMO FINAL DOS RESULTADOS');

% Tabela 1
disp('Tabela 1 (Var mu, P=1)');
if isempty(H1_table)
    disp('Nenhuma solução viável encontrada para o Teste 1.');
else
    T1 = array2table(H1_table, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'Max_Gap', 'Variaveis', 'Linhas_LMI'});
    disp(T1);
    writetable(T1, 'sistema_teste_d_Hinf_Tabela_H.csv');
end

% Tabela 2
disp('Tabela 2 (Var P, mu = 0, 1 e 8)');
if isempty(H2_table)
    disp('Nenhuma solução viável encontrada para o Teste 2.');
else
    T2 = array2table(H2_table, 'VariableNames', {'Grau_P', ...
        'Erro_mu0', 'MaxGap_mu0', 'V_mu0', 'L_mu0', ...
        'Erro_mu1', 'MaxGap_mu1', 'V_mu1', 'L_mu1', ...
        'Erro_mu8', 'MaxGap_mu8', 'V_mu8', 'L_mu8'});
    disp(T2);
    writetable(T2, 'sistema_teste_d_Hinf_Tabela_Hp.csv');
end

% Tabela 3
disp('Tabela 3 (Var mu, P=2)');
if isempty(H3_table)
    disp('Nenhuma solução viável encontrada para o Teste 3.');
else
    T3 = array2table(H3_table, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'Max_Gap', 'Variaveis', 'Linhas_LMI'});
    disp(T3);
    writetable(T3, 'sistema_teste_d_Hinf_Tabela_Hnovo.csv');
end

% Salva os gráficos
if N == 2
    if ~isempty(H1_table), figure(fig1); print('sistema_teste_d_Hinf_mu_0_a_10_P_1', '-depsc'); end
    if ~isempty(H2_table)
        figure(fig2); print('sistema_teste_d_Hinf_P_1_a_5_mu_0_1_8', '-depsc'); 
        figure(fig_bar); print('sistema_teste_d_Hinf_barras_P_1_a_5_mu_0_1_8', '-depsc');
    end
    if ~isempty(H3_table), figure(fig3); print('sistema_teste_d_Hinf_mu_pares_ate_14_P_2', '-depsc'); end
    fprintf('Gráficos e tabelas salvos no diretório atual.');
else
    if ~isempty(H2_table)
        % O gráfico de barras não depende de N=2, então podemos salvá-lo
        figure(fig_bar); print('sistema_teste_d_Hinf_barras_P_1_a_5_mu_0_1_8', '-depsc');
        fprintf('Gráfico de barras e tabelas salvos. Gráficos de linha ignorados (N=%d).', N);
    end
end
fprintf('Teste concluído.');