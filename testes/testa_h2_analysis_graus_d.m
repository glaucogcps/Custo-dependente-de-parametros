% function testa_h2_analysis()

close all
clear variables; 
load sistema_teste_d.mat;

% A norma H2 é infinita se D não for nula.
for k = 1:numel(D)
    D{k} = zeros(size(D{k}));
end
disp('Matrizes D foram zeradas para cálculo da norma H2.');
vetor_cor=['r--';'b--';'g--';'c--';'m--';'b-.';'g-.';'r-.';'c-.';'m-.'];
%  Teste 1: Variando o grau de Gamma (degGamma) com P fixo (deg=1) 
figure;
hold on;
title('Variação do grau (d) de \mu(\alpha) (P grau 1)');
ylabel('Norma H_2');
xlabel('\alpha');
H = [];

opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; 

disp('Iniciando Teste 1...');
for d = 0:8
    opt.deg = 1;        
    opt.degGamma = d;   
    
    out = h22_lmi_d_incerto_param(A, B, C, D, opt);
    
    if d == 0
        plot(out.alpha, out.realCosts, 'k', 'LineWidth', 1.5);
        plot(out.alpha, out.gcosts, 'b--', 'LineWidth', 1.5);
        
    else
        plot(out.alpha, out.gcosts, vetor_cor(d,:));
    end
    
    H = [H; d, norm(out.gcosts - out.realCosts), out.V, out.L];
end

disp('Tabela H (degGamma variando):');
disp(H);
legend('H_2', 'd=0', 'd=1', 'd=2', 'd=3', 'd=4', 'd=5', 'd=6', 'd=7', 'd=8', 'd=9', 'd=10', 'Location', 'southeast');
% Salva o gráfico 1
print('sistema_teste_d_H2_mu_0_a_10_P_1', '-depsc');

% Salva a Tabela H em CSV
T_H = array2table(H, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'V', 'L'});
writetable(T_H, 'sistema_teste_d_H2_Tabela_H.csv');


%  Teste 2: Variando o grau de P (deg) comparando Gamma fixo (0, 1 e 8) 
figure;
hold on;
title('Variação do grau (g) de P');
ylabel('Norma H_2');
xlabel('\alpha');
Hp = [];
it = 0;

disp('Iniciando Teste 2...');
for d = 1:5
    it = it + 1;
    
    % Caso 1: degGamma = 0
    opt.deg = d;
    opt.degGamma = 0;
    out1 = h22_lmi_d_incerto_param(A, B, C, D, opt);
    erro1 = norm(out1.gcosts - out1.realCosts);
    
    % Caso 2: degGamma = 1 
    opt.degGamma = 1;
    out2 = h22_lmi_d_incerto_param(A, B, C, D, opt);
    erro2 = norm(out2.gcosts - out2.realCosts);
    
    % Caso 3: degGamma = 8 (Alterado de 14 para 8)
    opt.degGamma = 8;
    out3 = h22_lmi_d_incerto_param(A, B, C, D, opt);
    erro3 = norm(out3.gcosts - out3.realCosts);
    
    if it == 1
        plot(out1.alpha, out1.realCosts, 'k', 'LineWidth', 1.5);
    end
    plot(out1.alpha, out1.gcosts, vetor_cor(d,:)); % Plota as curvas do Caso 1
    
    % Armazena métricas de todos os casos
    Hp = [Hp; d, erro1, out1.V, out1.L, erro2, out2.V, out2.L, erro3, out3.V, out3.L];
end

disp('Tabela Hp (degP variando, comparando degGamma 0, 1 e 8):');
disp(Hp);
legend('H_2', 'g=1', 'g=2', 'g=3', 'g=4', 'g=5', 'Location', 'southeast');
% Salva o gráfico 2
print('sistema_teste_d_H2_P_1_a_5_mu_0_1_8', '-depsc');

% Salva a Tabela Hp em CSV
T_Hp = array2table(Hp, 'VariableNames', {'Grau_P', 'Erro_mu0', 'V_mu0', 'L_mu0', 'Erro_mu1', 'V_mu1', 'L_mu1', 'Erro_mu8', 'V_mu8', 'L_mu8'});
writetable(T_Hp, 'sistema_teste_d_H2_Tabela_Hp.csv');

% =========================================================================
%  INSERÇÃO DO GRÁFICO DE BARRAS AGRUPADAS
% =========================================================================
grau = Hp(:,1);
diferenca_caso1 = Hp(:,2); % Erro com degGamma = 0
diferenca_caso2 = Hp(:,5); % Erro com degGamma = 1
diferenca_caso3 = Hp(:,8); % Erro com degGamma = 8

figure;
bar(grau, [diferenca_caso1, diferenca_caso2, diferenca_caso3]);
xlabel('Grau de P (g)');
ylabel('Norma da diferença (Custo Garantido - H_2)');
title('Comparação: Diferença vs Grau de P');
legend('\mu(\alpha) grau 0', '\mu(\alpha) grau 1', '\mu(\alpha) grau 8', 'Location', 'best');
grid on;
% Salva o gráfico de barras
print('sistema_teste_d_H2_barras_P_1_a_5_mu_0_1_8', '-depsc');
% =========================================================================


%  Teste 3: Variando Gamma (degGamma pares) com P fixo (deg=2) 
figure;
hold on;
title('Variação do grau (d) de \mu(\alpha) (P grau 2)');
ylabel('Norma H_2');
xlabel('\alpha');
Hnovo = [];
it = 0;

opt.deg = 2;

disp('Iniciando Teste 3...');
for d = 0:8       
    opt.degGamma = d;   
    
    out = h22_lmi_d_incerto_param(A, B, C, D, opt);
    
    if d == 0
        plot(out.alpha, out.realCosts, 'k', 'LineWidth', 1.5);
        plot(out.alpha, out.gcosts, 'b--', 'LineWidth', 1.5);
        
    else
        plot(out.alpha, out.gcosts, vetor_cor(d,:));
    end
    
    H = [H; d, norm(out.gcosts - out.realCosts), out.V, out.L];
end

disp('Tabela Hnovo (degGamma pares, degP=2):');
disp(Hnovo);
legend('H_2', 'd=1', 'd=2', 'd=3', 'd=4', 'd=5', 'd=6', 'd=7', 'd=8', 'Location', 'southeast');
% Salva o gráfico 3
print('sistema_teste_d_H2_mu_pares_ate_14_P_2', '-depsc');

% Salva a Tabela Hnovo em CSV
T_Hnovo = array2table(Hnovo, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'V', 'L'});
writetable(T_Hnovo, 'sistema_teste_d_H2_Tabela_Hnovo.csv');