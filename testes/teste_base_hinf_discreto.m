clear; clc;
load('polytop_hinf_d_finite_grid_vs_lmi_N2_n3_step0.01.mat');

% tab = [];
% for i=1:10
%     A = matrix_configs{i,1};
%     B = matrix_configs{i,2};
%     C = matrix_configs{i,3};
%     D = matrix_configs{i,4};
%     dados = matrix_configs{i, 5};
%     for g = 0:5
%         out = hinf_lmi_d_incerto(A,B,C,D, 'grau', g);
%         tab = [tab; i, g, out.hinf, dados.max_hinf_norm, abs(out.hinf-dados.max_hinf_norm), 100*abs(out.hinf-dados.max_hinf_norm)/dados.max_hinf_norm] 
%     end
% end

tab = [];
for i=5
    A = matrix_configs{i,1};
    B = matrix_configs{i,2};
    C = matrix_configs{i,3};
    D = matrix_configs{i,4};
    dados = matrix_configs{i, 5};
    for g = 1
        out = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 1);
        output = hinf_sistema_incerto_discrete(A, B, C, D, 'hinf', out.hinf, 'plot_graphic', false);
        %tab = [tab; i, g, out.hinf, dados.max_hinf_norm, abs(out.hinf-dados.max_hinf_norm), 100*abs(out.hinf-dados.max_hinf_norm)/dados.max_hinf_norm] 
    end
end