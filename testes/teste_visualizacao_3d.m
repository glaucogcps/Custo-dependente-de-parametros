
num_dimensoes = 3;
passo = 0.05; % Malha
alphas_3d = particao_simplex(num_dimensoes, passo);

alpha_matrix = cell2mat(alphas_3d); 

figure;
plot3(alpha_matrix(1,:), alpha_matrix(2,:), alpha_matrix(3,:), 'b.', 'MarkerSize', 8);
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('\alpha_3');
title('Politopo (N=3 Vértices)');
grid on;
axis equal;
view(3); % Define a vista 3D
rotate3d on; % Permite rotacionar o gráfico com o mouse