% function testa_hinf_analysis()

close all
clear variables;
load sistema_teste_d.mat;

%  Teste 1: Variando o grau de Gamma (degGamma) com P fixo (deg=1) 
figure;
hold on;
title('Variação do grau (d) de \mu(\alpha) (P grau 1)');
ylabel('Norma H_\infty');
xlabel('\alpha');
H = [];

% Configurações fixas iniciais
opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; % Usando a formulação relaxada (Finsler)
vetor_cor=['r--';'b--';'g--';'c--';'m--';'b-.';'g-.';'r-.';'c-.';'m-.'];
disp('Iniciando Teste 1...');
for d = 0:10
    % Atualiza opções para a nova função
    opt.deg = 1;        % Antigo degP
    opt.degGamma = d;   % Antigo degGamma
    
    out = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    
    if d == 0
        % Plota custo garantido (azul) e custo real (vermelho)
        plot(out.alpha, out.realCosts, 'k', 'LineWidth', 1.5);
        plot(out.alpha, out.gcosts, 'b--', 'LineWidth', 1.5);
        
    else
        plot(out.alpha, out.gcosts, vetor_cor(d,:));
    end
    
    % Armazena métricas: [grau, erro_norma, num_vars, num_linhas_LMI]
    H = [H; d, norm(out.gcosts - out.realCosts), out.V, out.L];
end

disp('Tabela H (degGamma variando):');
disp(H);
legend('H_\infty', 'd=0', 'd=1', 'd=2', 'd=3', 'd=4', 'd=5', 'd=6', 'd=7', 'd=8', 'd=9', 'd=10', 'Location', 'southeast');
% Salva o gráfico 1
print('sistema_teste_d_Hinf_mu_0_a_10_P_1', '-depsc');

% Salva a Tabela H em CSV
T_H = array2table(H, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'V', 'L'});
writetable(T_H, 'sistema_teste_d_Hinf_Tabela_H.csv');


%  Teste 2: Variando o grau de P (deg) comparando Gamma fixo (0, 1 e 8) 
figure;
hold on;
title('Variação do grau (g) de P');
ylabel('Norma H_\infty');
xlabel('\alpha');
Hp = [];
it = 0;

disp('Iniciando Teste 2...');
for d = 1:5
    it = it + 1;
    
    % Caso 1: degGamma = 0
    opt.deg = d;
    opt.degGamma = 0;
    out1 = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    erro1 = norm(out1.gcosts - out1.realCosts);
    
    % Caso 2: degGamma = 1 
    opt.degGamma = 1;
    out2 = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    erro2 = norm(out2.gcosts - out2.realCosts);
    
    % Caso 3: degGamma = 8 (Alterado de 14 para 8)
    opt.degGamma = 8;
    out3 = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    erro3 = norm(out3.gcosts - out3.realCosts);
    
    if it == 1
        plot(out1.alpha, out1.realCosts, 'k', 'LineWidth', 1.5);
    end
    plot(out1.alpha, out1.gcosts, vetor_cor(d,:));
    
    % Armazena métricas de todos os casos
    Hp = [Hp; d, erro1, out1.V, out1.L, erro2, out2.V, out2.L, erro3, out3.V, out3.L];
end

disp('Tabela Hp (degP variando, comparando degGamma 0, 1 e 8):');
disp(Hp);
legend('H_\infty', 'g=1', 'g=2', 'g=3', 'g=4', 'g=5', 'Location', 'southeast');
% Salva o gráfico 2
print('sistema_teste_d_Hinf_P_1_a_5_mu_0_1_8', '-depsc');

% Salva a Tabela Hp em CSV (Atualizado para mu8)
T_Hp = array2table(Hp, 'VariableNames', {'Grau_P', 'Erro_mu0', 'V_mu0', 'L_mu0', 'Erro_mu1', 'V_mu1', 'L_mu1', 'Erro_mu8', 'V_mu8', 'L_mu8'});
writetable(T_Hp, 'sistema_teste_d_Hinf_Tabela_Hp.csv');

% 
%  INSERÇÃO DO GRÁFICO DE BARRAS AGRUPADAS
% 
grau = Hp(:,1);
diferenca_caso1 = Hp(:,2); % Erro com degGamma = 0
diferenca_caso2 = Hp(:,5); % Erro com degGamma = 1
diferenca_caso3 = Hp(:,8); % Erro com degGamma = 8 (Coluna 8 da matriz Hp)

figure; % Abre uma nova janela especificamente para o gráfico de barras
bar(grau, [diferenca_caso1, diferenca_caso2, diferenca_caso3]);
xlabel('Grau de P (g)');
ylabel('Norma da diferença (Custo Garantido - H_\infty)');
title('Comparação H_\infty: Diferença vs Grau de P');
legend('\mu(\alpha) grau 0', '\mu(\alpha) grau 1', '\mu(\alpha) grau 8', 'Location', 'best');
grid on;
% Salva o gráfico de barras
print('sistema_teste_d_Hinf_barras_P_1_a_5_mu_0_1_8', '-depsc');
%


%  Teste 3: Variando Gamma (degGamma pares) com P fixo (deg=2) 
figure;
hold on;
title('Variação do grau (d) de \mu(\alpha) (P grau 2)');
ylabel('Norma H_\infty');
xlabel('\alpha');
Hnovo = [];
it = 0;

% Fixa grau de P em 2
opt.deg = 2;

disp('Iniciando Teste 3...');
for d = 1:8
    it = it + 1;
    
    % Atualiza grau de Gamma
    opt.degGamma = d;
    
    % Chamada da nova função
    out = hinf_lmi_d_incerto_param(A, B, C, D, opt);
    
    if it == 1
        plot(out.alpha, out.realCosts, 'k', 'LineWidth', 1.5);
    end
    plot(out.alpha, out.gcosts, vetor_cor(d,:));
    
    Hnovo = [Hnovo; d, norm(out.gcosts - out.realCosts), out.V, out.L];
end

disp('Tabela Hnovo (degGamma pares, degP=2):');
disp(Hnovo);
legend('H_\infty', 'd=1', 'd=2', 'd=3', 'd=4', 'd=5', 'd=6', 'd=7', 'd=8', 'Location', 'southeast');
% Salva o gráfico 3
print('sistema_teste_d_Hinf_mu_pares_ate_14_P_2', '-depsc');

% Salva a Tabela Hnovo em CSV
T_Hnovo = array2table(Hnovo, 'VariableNames', {'Grau_mu', 'Erro_Norma', 'V', 'L'});
writetable(T_Hnovo, 'sistema_glauco_Hinf_Tabela_Hnovo.csv');
