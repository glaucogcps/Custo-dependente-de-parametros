function P = gen_coefs(N, d)
% function P = gen_coefs(N, d)
%
% Gera uma matriz contendo os expoentes de todos os monômios de grau 'd'
% com 'N' variáveis.
%
% Esta função é utilizada para construir polinômios homogêneos em estruturas
% de otimização robusta (ROLMIP), onde a soma dos expoentes em cada linha
% deve ser igual a 'd'.
%
% Input:
%   N -> Número de variáveis (dimensão do simplex/número de vértices).
%   d -> Grau do polinômio.
%
% Output:
%   P -> Matriz de dimensão (n_monomios x N), onde cada linha representa
%        um monômio e cada coluna o expoente da respectiva variável.
%        O número de linhas é dado por n_monomios = (N+d-1)! / (d! * (N-1)!).
%
% Exemplo:
%   gen_coefs(3, 2) para variáveis x, y, z retorna:
%   2 0 0  (x^2)
%   1 1 0  (xy)
%   1 0 1  (xz)
%   0 2 0  (y^2)
%   0 1 1  (yz)
%   0 0 2  (z^2)

    if d == 0
        P = zeros(1, N);
        return;
    end

    if N == 1
        P = d;
        return;
    end

    % Algoritmo recursivo para gerar as partições
    P = [];
    for i = d:-1:0
        % Fixa o expoente da primeira variável como 'i'
        % e chama recursivamente para as N-1 variáveis restantes
        % com o grau restante (d-i).
        temp = gen_coefs(N-1, d-i);
        
        % Concatena o expoente fixo 'i' com as combinações retornadas
        block = [i * ones(size(temp, 1), 1), temp];
        P = [P; block];
    end
end