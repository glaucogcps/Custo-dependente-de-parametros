fig1 = openfig('teste_controle_H2_Ex3_N2_grafico.fig', 'reuse');
fig2 = openfig('teste_controle_Hinf_Ex3_N2_grafico.fig', 'reuse');

ax1_original = findobj(fig1, 'type', 'axes', '-depth', 1);
ax2_original = findobj(fig2, 'type', 'axes', '-depth', 1);

fig_nova = figure('Position', [100 100 1400 500]);

% ===== SUBPLOT 1 =====
ax1_novo = subplot(1, 2, 1);

% Copiar TODOS os filhos (linhas, patches, etc.)
copyobj(ax1_original.Children, ax1_novo);

% Copiar propriedades do axes
ax1_novo.XLim = ax1_original.XLim;
ax1_novo.YLim = ax1_original.YLim;
ax1_novo.XLabel.String = ax1_original.XLabel.String;
ax1_novo.YLabel.String = ax1_original.YLabel.String;
ax1_novo.Title.String = ax1_original.Title.String;
ax1_novo.XGrid = ax1_original.XGrid;
ax1_novo.YGrid = ax1_original.YGrid;

% Copiar legenda MANUALMENTE (se existir)
if ~isempty(ax1_original.Legend)
    leg_strings = ax1_original.Legend.String;  % Pega os textos da legenda
    legend(ax1_novo, leg_strings, 'Location', 'southwest', 'Orientation', 'vertical');
end

% ===== SUBPLOT 2 =====
ax2_novo = subplot(1, 2, 2);

% Copiar TODOS os filhos
copyobj(ax2_original.Children, ax2_novo);

% Copiar propriedades do axes
ax2_novo.XLim = ax2_original.XLim;
ax2_novo.YLim = ax2_original.YLim;
ax2_novo.XLabel.String = ax2_original.XLabel.String;
ax2_novo.YLabel.String = ax2_original.YLabel.String;
ax2_novo.Title.String = ax2_original.Title.String;
ax2_novo.XGrid = ax2_original.XGrid;
ax2_novo.YGrid = ax2_original.YGrid;

% Copiar legenda MANUALMENTE (se existir)
if ~isempty(ax2_original.Legend)
    leg_strings = ax2_original.Legend.String;  % Pega os textos da legenda
    legend(ax2_novo, leg_strings, 'Location', 'southwest', 'Orientation', 'vertical');
end

% Ajustar espaçamento entre subplots
% Reduz o espaço horizontal entre eles
subplot(1, 2, 1);
pos1 = get(gca, 'Position');
set(gca, 'Position', [pos1(1) pos1(2) pos1(3)*1 pos1(4)]);

subplot(1, 2, 2);
pos2 = get(gca, 'Position');
set(gca, 'Position', [pos2(1)-0.05 pos2(2) pos2(3)*1 pos2(4)]);


% Salvar como .fig
savefig(fig_nova, 'combinada.fig');

% Salvar como .eps (alta qualidade)
print(fig_nova, 'combinada.eps', '-depsc');

% Fechar figuras originais
close(fig1);
close(fig2);
