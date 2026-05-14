function out = h2_lmi_c_incerto_param(A, B, C, D, varargin)
% function out = h2_lmi_c_incerto_param(A, B, C, D, varargin)
%
% Avalia o custo H2 garantido de um SISTEMA LINEAR CONTÍNUO sujeito a 
% incertezas politópicas, utilizando funções de Lyapunov dependentes de parâmetros.
%
% Entradas:
%   A, B, C, D - Cell arrays contendo as matrizes dos vértices do politopo.
%                Ex: A = {A1, A2, ..., AN}, onde N é o número de vértices.
%                Nota: Para norma H2 contínua finita, a matriz D deve ser nula.
%   varargin   - (Opcional) Struct com opções de configuração:
%                .solver   : Solver a ser utilizado (default: 'mosek')
%                .deg      : Grau do polinômio da matriz de Lyapunov P (default: 1)
%                .degGamma : Grau do polinômio da variável de desempenho (default: 0)
%                .op       : Tipo de otimização (0 = integral/média, 1 = pico) (default: 0)
%                .varFolga : Uso de variáveis de folga/Finsler (1 = sim, 0 = não) (default: 1)
%                .verbose  : Nível de detalhe do solver (default: 0)
%
% Saída:
%   out - Struct contendo os resultados da otimização e validação:
%         .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%         .wc        : Pior custo H2 garantido encontrado no grid
%         .gcosts    : Vetor com os custos garantidos avaliados no grid
%         .realCosts : Vetor com os custos reais (norma 2) avaliados no grid
%         .alpha     : Matriz com os pontos do simplex testados (cada linha é um ponto)
%         .L         : Número de linhas das LMIs
%         .V         : Número de variáveis escalares
%         .cpusec_m  : Tempo de montagem das LMIs
%         .cpusec_s  : Tempo de resolução do solver
%
% Date: 22/03/2026

%   Validação  
if nargin < 4, error('Faltam argumentos'); end
num_vertices = numel(A); 
n = size(A{1}, 1); m = size(B{1}, 2); p = size(C{1}, 1);

%   Opções  
options = struct('solver', 'mosek', 'deg', 1, 'degGamma', 0, 'op', 0, 'varFolga', 1, 'verbose', 0);
if nargin > 4 && isstruct(varargin{1})
    f = fieldnames(varargin{1});
    for i = 1:length(f), options.(f{i}) = varargin{1}.(f{i}); end
end

%   Variáveis ROLMIP  
A_rol = rolmipvar(A, 'A', num_vertices, 1);
B_rol = rolmipvar(B, 'B', num_vertices, 1);
C_rol = rolmipvar(C, 'C', num_vertices, 1);

P = rolmipvar(n, n, 'P', 'symmetric', num_vertices, options.degP);

% Configuração da variável de desempenho (rho/gamma)
if options.op == 1
    % Caso op=1: rho é uma variável escalar constante (pico)
    rho = sdpvar(1, 1);
    rho_alpha = rolmipvar(rho, 'rho', num_vertices, 0); % Converte para rolmipvar grau 0
    obj = rho;
else        
    % Caso op=0: rho é um polinômio rho(alpha)
    powers = gen_coefs(num_vertices, options.degGamma);
    obj = [];
    rhos = cell(1, size(powers, 1));
    
    for i = 1:size(powers, 1)
        rhos{i} = sdpvar(1, 1);
        % Coeficiente binomial para ponderação na função objetivo (integral/média)
        c_val = factorial(options.degGamma) / prod(factorial(powers(i, :)));  
        if isempty(obj)
            obj = rhos{i} * (1/c_val);
        else
            obj = rhos{i} * (1/c_val) + obj;
        end
    end    
    rho_alpha = rolmipvar(rhos, 'rho', num_vertices, options.degGamma);
end

%   LMIs  
if options.varFolga == 0  
    LMIs = [trace(B_rol'*P*B_rol) <= rho_alpha, P >= 0, A_rol'*P + P*A_rol + C_rol'*C_rol <=0];
else
    LMIs = [trace(B_rol'*P*B_rol) <= rho_alpha];
    
    %   Formulação Relaxada (Variáveis de Folga / Finsler)  
    % Introduz multiplicador X para desacoplar as matrizes do sistema de P
    
    X = rolmipvar(2*n + p, n + p, 'X', 'full', num_vertices, options.deg);
    Z = zeros(n);
    
    % Matriz Q contendo P e rho
    Q = blkdiag([Z, P; P, Z], eye(p));
    
    % Matriz Bcal contendo a dinâmica do sistema
    % Estrutura típica para impor a dinâmica via multiplicadores
    Bcal = [eye(n), -A_rol, zeros(n, p);
            zeros(p, n), -C_rol, eye(p)];
            
    % LMI: Q + X*Bcal + Bcal'*X' <= 0
    LMIs = [LMIs, Q + X * Bcal + Bcal' * X' <= 0, P >= 0];
end

%   Solução  
info = lmiinfo(LMIs);
out.L = sum(info.sdp(:, 1));
out.V = size(getvariables(LMIs), 2);

out.cpusec_m = clock;
sol = solvesdp(LMIs, obj, sdpsettings('verbose', options.verbose, 'solver', options.solver));
out.cpusec_s = sol.solvertime;

if min(checkset(LMIs)) > -1e-6
    out.feas = 1;
    
    % Calcula pior caso (Grid) generalizado para N vértices
    max_val = -inf;
    out.gcosts = [];
    out.realCosts = [];
    
    % Define o passo da malha (pode ser parametrizado no varargin futuramente)
    passo_grid = 0.05; 
    
    % Gera os pontos do simplex unitário usando a nova função
    pontos_alpha = particao_simplex(num_vertices, passo_grid);
    
    % Inicializa a matriz para salvar os pontos testados (útil para plotagem)
    out.alpha = zeros(length(pontos_alpha), num_vertices);
    
    % Preserva originais para validação
    Ao = A; Bo = B; Co = C;
    
    for k = 1:length(pontos_alpha)
        alpha_vec = pontos_alpha{k}; % Vetor coluna com as coordenadas baricêntricas
        out.alpha(k, :) = alpha_vec'; % Salva como linha na saída
        
        % 1. Custo Garantido
        % evalpar do ROLMIP espera um cell array com o vetor de coordenadas (em linha)
        W_val = double(evalpar(rho_alpha, {alpha_vec'}));
        v = sqrt(trace(W_val));
        out.gcosts = [out.gcosts; v];
        if v > max_val, max_val = v; end
        
        % 2. Custo Real
        % Inicializa matrizes vazias para a combinação convexa
        Aa = zeros(n, n);
        Ba = zeros(n, m);
        Ca = zeros(p, n);
        
        % Combinação convexa generalizada para N vértices
        for i = 1:num_vertices
            Aa = Aa + alpha_vec(i) * Ao{i};
            Ba = Ba + alpha_vec(i) * Bo{i};
            Ca = Ca + alpha_vec(i) * Co{i};
        end
        
        % Monta o sistema contínuo (D = 0 para norma H2 finita)
        sys = ss(Aa, Ba, Ca, zeros(p,m)); 
        try
            out.realCosts = [out.realCosts; norm(sys, 2)];
        catch
            out.realCosts = [out.realCosts; NaN];
        end
    end
    out.wc = max_val;
else
    out.feas = 0;
    out.rho_alpha = [];
    out.wc = Inf;
    out.gcosts = [];
    out.realCosts = [];
    out.alpha = [];
end
end