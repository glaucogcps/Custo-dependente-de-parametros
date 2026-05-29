# Dependências
import pandas            as pd
import matplotlib as mpl
mpl.use("pgf")
mpl.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,
    "font.family": "serif",
    "pgf.rcfonts": False,
    "font.size": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pgf.preamble": r"""
        \usepackage{amsmath}
        \usepackage{amssymb}
        \usepackage{mathtools}
        \usepackage{dsfont}
        \providecommand{\mathdefault}[1]{#1}
    """
})
import matplotlib.pyplot as plt

def plot_from_csv(file,
                  xlabel,
                  ylabel,
                  out,
                  xlim,
                  pmin=None,
                  pmax=None):

    # Lê o CSV
    df = pd.read_csv(file)

    # Ordena para garantir curvas corretas
    df = df.sort_values(by=["Grau_rho", "Alpha"])

    # Cria figura
    fig, ax = plt.subplots(figsize=(3.4, 2.2))

    # Itera sobre cada grau
    for grau, grupo in df.groupby("Grau_rho"):

        alpha = grupo["Alpha"].values

        # --------------------------------------------
        # Transformação opcional do eixo x
        # --------------------------------------------
        if pmin is not None and pmax is not None:
            x = pmax * alpha + pmin * (1 - alpha)
        else:
            x = alpha

        custo_g = grupo["Custo_Garantido"].values
        custo_r = grupo["Custo_Real"].values

        # Linha contínua -> custo real
        line_real, = ax.plot(
            x,
            custo_r,
            linewidth=1.5,
            label=fr"$d={grau}$"
        )

        # Usa mesma cor da linha anterior
        cor = line_real.get_color()

        # Linha tracejada -> custo garantido
        ax.plot(
            x,
            custo_g,
            "--",
            color=cor,
            linewidth=1.5
        )

    # Configurações do gráfico
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.grid(
        True,
        alpha=0.25,
        linestyle="--",
        linewidth=0.5
    )

    ax.set_xlim(xlim)

    # Ajusta espaço interno para acomodar legenda
    fig.subplots_adjust(right=0.73)

    # Legenda externa à direita
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),

        frameon=False,
        fancybox=False,

        handlelength=0.8,
        handletextpad=0.4,
        borderaxespad=0.2,
        labelspacing=0.3
    )

    # Layout final
    plt.tight_layout()

    # Salva em PGF para LaTeX
    plt.savefig(out, bbox_inches="tight")
    plt.close
