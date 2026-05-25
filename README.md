# Custo dependente de parâmetros

Este repositório trata de uma pesquisa de IC focada em controle robusto $H_2$ e $H_\infty$ de sistemas incertos politópicos. O foco do projeto é tratar a variável de otimização custo garantido como uma $\text{rolmipvar}$, ou seja, o custo garantido também depende da variável de incerteza.
Testes demonstraram que essa abordagem permite uma melhor aproximação entre as curvas de custo garantido e da norma verdadeira ao longo do politopo.

## Estrutura do repositório

O repositório apresenta duas pastas principais: $\text{funcoes}$ e $\text{testes}$. A primeira é dividida em outras subpastas, cada qual com uma responsabilidade específica:
$\text{LMIs}$ trata de condições clássicas de análise para sistemas precisamente conhecidos e com incerteza politópica (com cada código utilizando uma formulação LMI diferente), 
$\text{controle}$ trata de funções de controle de SLIT com incerteza politópica utilizando formulações de LMI, $\text{MJLS}$ possui funções de análise para sistemas
MJLS precisamente conhecidos e com incerteza politópica na matriz de probabilidades e controle de MJLS com incerteza na matriz de probabilidades, e 
$\text{auxiliares}$ possui diversas funções secundárias da apoio (das quais vale destacar "particao_simplex.m" que gera e salva os pontos de um grid uniforme
ao longo do simplex unitário. Todas as funções são testadas na pasta de $\text{testes}$.

## Modo de uso
Uma vez baixado o repositório, adicione a pasta $\text{funcoes}$ ao PATH do MATLAB. Em seguida, ao abrir a pasta $\text{testes}$, escolha uma funcionalidade 
específica para testar. Coloque a pasta "alfas_simplex" dentro do repositório no qual o script de teste está rodando para acelerar o código (isso será consertado posteriormente).

## Referências (resumido)
[1] K. Zhou, J. C. Doyle, and K. Glover, Robust and Optimal Control.
Upper Saddle River, NJ, USA: Prentice Hall, 1996. </br>
[2] S. Boyd, L. El Ghaoui, E. Feron, and V. Balakrishnan, Linear Matrix
Inequalities in System and Control Theory. Philadelphia, PA: SIAM
Studies in Applied Mathematics, 1994.</br>
[3] MOSEK ApS, The MOSEK optimization software, 2015,
http://www.mosek.com.</br>
[4] O. L. V. Costa, M. D. Fragoso, and M. G. Todorov, Continuous-Time
Markov Jump Linear Systems. Berlin, Germany: Springer-Verlag, 2013.</br>
[5] O. L. V. Costa, M. D. Fragoso, and R. P. Marques, Discrete-Time
Markovian Jump Linear Systems. New York, NY, USA: Springer-
Verlag, 2005.</br>
[6] R. C. L. F. Oliveira, A. N. Vargas, J. B. R. do Val, and P. L. D. Peres,
“Robust stability, $H_2$ analysis and stabilisation of discrete-time Markov
jump linear systems with uncertain probability matrix,” Int. J. Control,
vol. 82, no. 3, pp. 470–481, Mar. 2009.</br>
[7] ——, “ $H_2$ guaranteed cost computation by means of parameter dependent
Lyapunov functions,” Int. J. Syst. Sci., vol. 35, no. 5, pp. 305–315,
Apr. 2004.</br>
[8] M. C. de Oliveira, J. C. Geromel, and L. Hsu, “LMI characterization of
structural and robust stability: The discrete-time case,” Lin. Alg. Appl.,
vol. 296, no. 1–3, pp. 27–38, Jun. 1999.</br>
[9] J. C. Geromel and R. H. Korogui, “Analysis and synthesis of robust
control systems using linear parameter dependent Lyapunov functions,”
IEEE Trans. Autom. Control, vol. 51, no. 12, pp. 1984–1989, Dec. 2006.</br>
[10] C. F. Morais, M. F. Braga, R. C. L. F. Oliveira, and P. L. D. Peres, “$H_2$
control of discrete-time Markov jump linear systems with uncertain transition
probability matrix: Improved linear matrix inequality relaxations
and multi-simplex modeling,” IET Control Theory & Appl., vol. 7, pp.
1665–1674, Aug. 2013.</br>
[11] R. C. L. F. Oliveira and P. L. D. Peres, “A convex optimization procedure
to compute H2 and H∞ norms for uncertain linear systems in polytopic
domains,” Optim. Control Appl. Meth., vol. 29, no. 4, pp. 295–312,
July/August 2008.</br>
[12] M. F. Braga, C. F. Morais, R. C. L. F. Oliveira, and P. L. D. Peres,
“Realimentação de estados robusta $H_2$ e $H_\infty$ para sistemas lineares
discretos no tempo por meio de LMIs com parâmetros escalares,” in XIX
CBA, Campina Grande, PB, Brasil, Setembro 2012, pp. 1672–1679.</br>
[13] C. F. Morais, M. F. Braga, R. C. L. F. Oliveira, and P. L. D. Peres, “Robust
$H_\infty$ state feedback control for MJLS with uncertain probabilities,”
Automatica, vol. 52, no. 2, pp. 317–431, Feb. 2015.</br>
[14] R. C. L. F. Oliveira and P. L. D. Peres, “Parameter-dependent LMIs in robust
analysis: Characterization of homogeneous polynomially parameterdependent
solutions via LMI relaxations,” IEEE Trans. Autom. Control,
vol. 52, no. 7, pp. 1334–1340, Jul. 2007.</br>
[15] C. M. Agulhari, A. Felipe, R. C. L. F. Oliveira, and P. L. D. Peres,
“Algorithm 998: The Robust LMI Parser — A toolbox to construct LMI
conditions for uncertain systems,” ACM Trans. Math. Softw., vol. 45,
no. 3, pp. 36:1–36:25, Aug. 2019, http://rolmip.github.io.</br>
[16] J. Löfberg, “YALMIP: A toolbox for modeling and optimization in MATLAB,”
in Proc. 2004 IEEE Int. Symp. on Comput. Aided Control Syst.
Des., Taipei, Taiwan, Sep. 2004, pp. 284–289, http://yalmip.github.io.
