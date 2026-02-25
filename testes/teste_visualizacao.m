num_dimensoes = 2;
passo = 0.01; % Malha fina
alphas_2d = particao_simplex(num_dimensoes, passo);

% Converter o cell array de vetores linha para uma matriz
alpha_matrix = cell2mat(alphas_2d); 

figure;
% Plotar as coordenadas X (posições ímpares) vs as coordenadas Y (posições pares)
plot(alpha_matrix(1,:), alpha_matrix(2,:), 'b.', 'MarkerSize', 8); % <--- CORREÇÃO AQUI
xlabel('\alpha_1');
ylabel('\alpha_2');
title('Politopo (N=2 Vértices)');
grid on;
axis equal; % Garante que os eixos tenham a mesma escala para uma visualização correta