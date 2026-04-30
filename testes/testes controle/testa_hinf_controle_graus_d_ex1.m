clear; clc; close all;

fprintf(' Teste de Síntese H-infinito - Exemplo 1 (de Oliveira et al.) \n');

%% 1. Definição do Sistema Precisamente Conhecido (N = 1)

nx = 4; % Número de estados (A é 4x4)
nu = 2; % Número de entradas de controle (B2 é 4x2)
nw = 3; % Número de entradas de distúrbio (B1 é 4x3)
nz = 3; % Número de saídas de desempenho (C é 3x4)
N  = 1; % Sistema precisamente conhecido (sem incerteza)

fprintf('Parâmetros: nx=%d, nu=%d, nw=%d, nz=%d, Vértices N=%d.\n', nx, nu, nw, nz, N);

A = cell(1, N); B = cell(1, N); E = cell(1, N);
Cz = cell(1, N); Dz = cell(1, N); Ez = cell(1, N);

% Inserindo as matrizes exatas do artigo

A{1} = [ 0.8189   0.0863   0.0900   0.0813;
         0.2524   1.0033   0.0313   0.2004;
        -0.0545   0.0102   0.7901  -0.2580;
        -0.1918  -0.1034   0.1602   0.8604];

B{1} = [ 0.0045   0.0044;
         0.1001   0.0100;
         0.0003  -0.0136;
        -0.0051   0.0936];

E{1} = [ 0.0953   0   0;
         0.0145   0   0;
         0.0862   0   0;
        -0.0011   0   0];

Cz{1} = [ 1   0  -1   0;
          0   0   0   0;
          0   0   0   0];

Dz{1} = [ 0   0;
          1   0;
          0   1];

% Matriz D11 (Ez) não fornecida no artigo, assumida como zero
Ez{1} = zeros(nz, nw); 

disp('Matrizes do Exemplo 1 carregadas com sucesso.');

fprintf('\n');

%% 2. Configuração das Opções do Solver

opt = struct();
opt.solver = 'mosek'; 
opt.verbose = 0;
opt.varFolga = 1; % Usa a formulação relaxada
opt.op = 1;       % Otimiza o pico do custo (rho constante) - Melhor para N=1

%% 3. Teste de Síntese

disp('Iniciando Teste de Síntese para o sistema nominal...');
% Como N=1, os graus polinomiais não têm efeito, podemos deixar como 0
opt.deg = 0;        
opt.degGamma = 0;   
% Chama a função de controle
out = hinf_lmi_d_incerto_controle(A, B, E, Cz, Dz, Ez, opt);
H_table = [];

if out.feas == 1
    fprintf('\n LMI Viável! (Tempo: %.2f s)\n', out.cpusec_s);
    % Como N=1, pegamos o primeiro (e único) elemento dos custos
    custo_garantido = out.gcosts(1);
    custo_real = out.realCosts(1);

    fprintf(' Custo H-infinito Garantido (LMI): %.4f\n', custo_garantido);
    fprintf(' Custo H-infinito Real (Malha Fechada): %.4f\n', custo_real);
    
    % Verifica se bate com o artigo
    if abs(custo_garantido - 1.248) < 0.01
        fprintf(' SUCESSO! O valor bateu com o g = 1.248 do artigo!\n');
    else
        fprintf(' O valor diferiu do artigo. Verifique a formulação da LMI.\n');
    end
    
    % Armazena dados para a tabela
    erro_norma = abs(custo_garantido - custo_real);
    max_gap = erro_norma;
    min_gap = erro_norma;
    norma_melhor_caso = custo_real;
    norma_pior_caso = custo_real;
    H_table = [opt.degGamma, erro_norma, max_gap, min_gap, norma_melhor_caso, norma_pior_caso, out.V, out.L];
else
    fprintf('\n LMI Inviável.\n');
end

%% 4. Exibição dos Resultados
disp(' Tabela de Resultados:');
if isempty(H_table)
    warning('Nenhuma solução viável foi encontrada.');
else
    % Cria a tabela com as novas colunas
    T_H = array2table(H_table, 'VariableNames', {'Grau_gamma', 'Erro_Norma', 'Max_Gap', 'Min_Gap', 'Norma_Melhor_Caso', 'Norma_Pior_Caso', 'Variaveis', 'Linhas_LMI'});
    disp(T_H);
    % Salva os resultados da tabela
    writetable(T_H, 'teste_controle_Hinf_Ex1_Tabela.csv');
    fprintf('Tabela salva como teste_controle_Hinf_Ex1_Tabela.csv\n');
end

% Salva o gráfico apenas se alguma figura estiver aberta
figHandles = findobj('Type', 'Figure');
if ~isempty(figHandles)
    savefig('teste_controle_Hinf_Ex1_grafico.fig');
    fprintf('Gráfico salvo como teste_controle_Hinf_Ex1_grafico.fig\n');
end

fprintf('Teste concluído.\n');