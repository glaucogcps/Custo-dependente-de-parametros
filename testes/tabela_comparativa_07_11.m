clear; clc; close all;
load('polytop_hinf_d_finite_grid_vs_lmi_N2_n3_step0.01.mat');
tab = [];
for i=5
    Ai = matrix_configs{i,1};
    Bi = matrix_configs{i,2};
    Ci = matrix_configs{i,3};
    Di = matrix_configs{i,4};
    prec = 100000;
    for j=1:2
        A{j} = double(int32(prec*Ai{j}))/prec;
        B{j} = double(int32(prec*Bi{j}))/prec;
        C{j} = double(int32(prec*Ci{j}))/prec;
        D{j} = double(int32(prec*Di{j}))/prec;
    end

    dados = matrix_configs{i, 5};
    output = hinf_sistema_incerto_discrete(A, B, C, D,  'plot_graphic', false);
    for g = 1:10
        out = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 0);
        h0 = out.hinf;
        % norma
        outmu = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 1);
        outmu2 = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 2);
        outmu3 = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 3);
        outmu8 = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 4);
        % soma coeffs
%         outmu = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 1);
%         outmu2 = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 2);
%         outmu3 = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 3);
%         outmu8 = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 4);
        H = [];
        for alpha=0:0.001:1
            if outmu.feas == 0
                h1 = 0;
            else
                h1 = sqrt(evalpar(outmu.hinf, {[alpha 1-alpha]}));
            end
            if outmu2.feas == 0
                h2 = 0;
            else
                h2 = sqrt(evalpar(outmu2.hinf, {[alpha 1-alpha]}));
            end
            if outmu3.feas == 0
                h3 = 0;
            else
                h3 = sqrt(evalpar(outmu3.hinf, {[alpha 1-alpha]}));
            end
            if outmu8.feas == 0
                h8 = 0;
            else
                h8 = sqrt(evalpar(outmu8.hinf, {[alpha 1-alpha]}));
            end
            Aa = alpha * A{1} + (1-alpha) * A{2};
            Ba = alpha * B{1} + (1-alpha) * B{2};
            Ca = alpha * C{1} + (1-alpha) * C{2};
            Da = alpha * D{1} + (1-alpha) * D{2};
            sys = ss(Aa, Ba, Ca, Da, -1);
            norma = norm(sys, inf);
            H = [H; alpha output.max_h_inf_norm h0 h1 h2 h3 h8 norma];
        end
         figure;
         hold on;
%         plot(H(:,1), H(:,2));
         plot(H(:,1), H(:,3));
         plot(H(:,1), H(:,4));
%         plot(H(:,1), H(:,5));
%         plot(H(:,1), H(:,6));
%         plot(H(:,1), H(:,7));
         plot(H(:,1), H(:,8));
%        % plot(output.tab(:,1), output.tab(:,3), '--');
         xlabel('\alpha');
         ylabel('H_\infty');
         %legend('norma', 'mu0', 'mu1', 'mu2', 'mu3', 'mu8', 'hinf');
         legend('\mu', '\mu(\alpha)', 'H_\infty');
         title(sprintf('Grau de P = %d', g));
%         
%          figure;
%          hold on;
% %         
%          plot(H(:,1), H(:,3) - H(:,8));
%          plot(H(:,1), H(:,4) - H(:,8));
%         plot(H(:,1), H(:,5) - H(:,8));
%         plot(H(:,1), H(:,6) - H(:,8));
%         plot(H(:,1), H(:,7) - H(:,8));
%        
%        % plot(output.tab(:,1), output.tab(:,3), '--');
%         xlabel('\alpha');
%         ylabel('H_\infty');
%         legend( '\mu - H_\infty', '\mu(\alpha) - H_\infty');
%         title(strcat('diff', num2str(g)));
        
        disp('g = ');
        disp(g);
        disp('norma da diferença grau 0');
        disp(norm( H(:,3) - H(:,8)));
        disp('num variaveis mu');
        disp(out.V)
        
        disp('norma da diferença grau 1');
        disp(norm( H(:,4) - H(:,8)));
        tab = [tab; g, norm(H(:,3) - H(:,8)), out.V, out.L, norm(H(:,4) - H(:,8)), outmu.V, outmu.L];
        
    end
    % Dados da tabela
    grau = tab(2:end,1);  % Grau de P

    % Diferenças para cada caso
    diferenca_caso1 = tab(2:end,2);
    diferenca_caso2 = tab(2:end,5);

    % Gráfico de barras agrupadas
    figure;
    bar(grau, [diferenca_caso1, diferenca_caso2]);
    xlabel('Grau de P (g)');
    ylabel('Norma da diferença');
    title('Comparação: Diferença vs Grau de P');
    legend('Caso 1', 'Caso 2');
    grid on;
end
