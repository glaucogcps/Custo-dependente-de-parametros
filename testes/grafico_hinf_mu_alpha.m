clear; clc; close all;
load('polytop_hinf_d_finite_grid_vs_lmi_N2_n3_step0.01.mat');

for i=5
    A = matrix_configs{i,1};
    B = matrix_configs{i,2};
    C = matrix_configs{i,3};
    D = matrix_configs{i,4};
    dados = matrix_configs{i, 5};
    output = hinf_sistema_incerto_discrete(A, B, C, D,  'plot_graphic', false);
    for g = 1:10
        out = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 0);
        h0 = out.hinf;
%         outmu = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 1);
%         outmu2 = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 2);
%         outmu3 = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 3);
%         outmu8 = hinf_lmi_d_incerto_schur(A,B,C,D, 'grau', g, 'grau_mu', 4);
        outmu = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 1);
        outmu2 = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 2);
        outmu3 = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 3);
        outmu8 = hinf_lmi_d_incerto(A,B,C,D, 'grau', g, 'grau_mu', 4);
        H = [];
        for alpha=0:0.01:1
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
        plot(H(:,1), H(:,2));
        plot(H(:,1), H(:,3));
        plot(H(:,1), H(:,4));
        plot(H(:,1), H(:,5));
        plot(H(:,1), H(:,6));
        plot(H(:,1), H(:,7));
        plot(H(:,1), H(:,8));
       % plot(output.tab(:,1), output.tab(:,3), '--');
        xlabel('alpha');
        ylabel('Hinf');
        legend('norma', 'mu0', 'mu1', 'mu2', 'mu3', 'mu8', 'hinf');
        title(num2str(g));
        
%          figure;
%         hold on;
%         
%         plot(H(:,1), H(:,3) - H(:,8));
%         plot(H(:,1), H(:,4) - H(:,8));
%         plot(H(:,1), H(:,5) - H(:,8));
%         plot(H(:,1), H(:,6) - H(:,8));
%         plot(H(:,1), H(:,7) - H(:,8));
       
       % plot(output.tab(:,1), output.tab(:,3), '--');
        xlabel('\alpha');
        ylabel('H_\infty');
        legend( 'mu0', 'mu1', 'mu2', 'mu3', 'mu8');
        title(strcat('diff', num2str(g)));
        
        disp('g = ');
        disp(g);
        disp('norma da diferença grau 0');
        disp(norm( H(:,3) - H(:,8)));
        disp('norma da diferença grau 1');
        disp(norm( H(:,4) - H(:,8)));
        disp('norma da diferença grau 2');
        disp(norm( H(:,5) - H(:,8)));
        disp('norma da diferença grau 3');
        disp(norm( H(:,6) - H(:,8)));
        disp('norma da diferença grau 8');
        disp(norm( H(:,7) - H(:,8)));
    end
end
