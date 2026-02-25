% SCRIPT DE GERAÇÃO FINAL (Com hinf_lmi_c_incerto_param)
clear; clc; close all;

fprintf('--- Geração de Politopos H-Infinito (Contínuo - Relaxado) ---\n');

% Configurações
n = 3; m = 2; p = 2; N = 2;
target_list_size = 10; 
simplex_step = 0.05; 

matrix_configs = {}; 
current_list_size = 0;
attempts = 0;

fprintf('Buscando %d sistemas...\n', target_list_size);
fprintf('%-10s %-15s %-15s %-15s\n', 'Tentativa', 'Status', 'Max Real Eig', 'Gap (%)');
fprintf('%s\n', repmat('-', 1, 60));

while current_list_size < target_list_size
    attempts = attempts + 1;
    
    % 1. Gera Sistema Nominal Estável
    A_nom = randn(n);
    max_real = max(real(eig(A_nom)));
    % Força estabilidade com margem de segurança (-0.5)
    if max_real > -0.5
        A_nom = A_nom - (max_real + 0.5 + rand()) * eye(n);
    end
    
    B_nom = randn(n, m);
    C_nom = randn(p, n);
    D_nom = randn(p, m);
    
    % 2. Gera Vértices por Perturbação (20% de variação)
    % Isso garante que os vértices sejam "parecidos" o suficiente para serem estáveis
    perturbation = 0.2; 
    
    A_cell = cell(1, N); B_cell = cell(1, N); C_cell = cell(1, N); D_cell = cell(1, N);
    stable_poly = true;
    
    for i = 1:N
        A_cell{i} = A_nom + randn(n)*perturbation;
        B_cell{i} = B_nom + randn(n, m)*perturbation;
        C_cell{i} = C_nom + randn(p, n)*perturbation;
        D_cell{i} = D_nom + randn(p, m)*perturbation;
        
        % Verifica estabilidade de cada vértice
        if max(real(eig(A_cell{i}))) >= -1e-4
            stable_poly = false; break;
        end
    end
    
    if ~stable_poly, continue; end

    % 3. Cálculo H-Infinito (LMI Relaxada)
    try
        % Usa a NOVA função que você forneceu
        % varFolga=1 (Relaxada) é crucial para encontrar soluções
        opt.varFolga = 1; 
        opt.deg = 1;      % Grau 1 para P(alpha)
        opt.verbose = 0;
        
        out_lmi = hinf_lmi_c_incerto_param(A_cell, B_cell, C_cell, D_cell, opt);
    catch
        out_lmi.feas = 0;
    end
    
    if out_lmi.feas == 0
        if mod(attempts, 50) == 0, fprintf('%d... (LMI inviável)\n', attempts); end
        continue;
    end
    
    % 4. Validação (Grid)
    out_grid = valor_singular_sistema_incerto_continuous(A_cell, B_cell, C_cell, D_cell, ...
        'simplex_step', simplex_step, 'plot_graphic', false);
    
    % 5. Critérios de Aceitação
    % Gap positivo e não trivial
    gap = (out_lmi.wc - out_grid.max_h_inf_norm) / out_grid.max_h_inf_norm;
    
    if gap > 0.001 
        current_list_size = current_list_size + 1;
        
        matrix_configs{current_list_size, 1} = A_cell;
        matrix_configs{current_list_size, 2} = B_cell;
        matrix_configs{current_list_size, 3} = C_cell;
        matrix_configs{current_list_size, 4} = D_cell;
        matrix_configs{current_list_size, 5} = struct(...
            'max_hinf_norm', out_grid.max_h_inf_norm, ...
            'alpha_star', out_grid.alpha_star, ...
            'custo_garantido', out_lmi.wc); % Note: out_lmi.wc na nova função
        
        fprintf('%-10d %-15s %-15.4f %-15.2f%%\n', attempts, 'SALVO!', max(real(eig(A_nom))), gap*100);
    end
end

filename = 'polytop_hinf_finite_grid_vs_lmi_N2_n3_step0.01_cont.mat';
save(filename, 'matrix_configs', 'N', 'n', 'm', 'p', 'simplex_step');
fprintf('Concluído. Arquivo salvo: %s\n', filename);