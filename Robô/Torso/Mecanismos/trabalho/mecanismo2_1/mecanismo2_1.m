% função x = 0 y = (l1 + l2)/2
%resultado conhecido = círculo
%escolha dos valores de forma a evitar singularidades
% w5 = -1 w6 = 1
%%Simulador de mecanismos planos
clear
%%
%Construção do modelo baseado em variáveis simbólicas
%Definição das variáveis
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4 x5 y5 phi5 x6 y6 phi6 rhof alphaf rhoe alphae real %posição
syms x1p y1p phi1p x2p y2p phi2p x3p y3p phi3p x4p y4p phi4p x5p y5p phi5p x6p y6p phi6p rhofp alphafp rhoep alphaep real %velocidade
syms l1 l2 %constantes do modelo
syms t %tempo

%Construção do vetor de coordenadas generalizadas e a sua derivada
q = [x1;y1;phi1;x2;y2;phi2;x3;y3;phi3;x4;y4;phi4;x5;y5;phi5;x6;y6;phi6;rhof;alphaf;rhoe;alphae];

qp = [x1p;y1p;phi1p;x2p;y2p;phi2p;x3p;y3p;phi3p;x4p;y4p;phi4p;x5p;y5p;phi5p;x6p;y6p;phi6p;rhofp;alphafp;rhoep;alphaep];

%Construção do vetor de restrições cinemáticas e diretoras
PhiK = [x1;
        y1;
        x2;
        y2;
        x5 - l1/2;
        y5;
        x6 + l1/2;
        y6;
        x1 - x3 + l1*cos(phi1);
        y1 - y3 + l1*sin(phi1);
        x2 - x4 - l1*cos(phi2);
        y2 - y4 - l1*sin(phi2);
        x3 - x4 - l2*sin(phi3) + l2*sin(phi4);
        y3 - y4 + l2*cos(phi3) - l2*cos(phi4);
        x5 - x3 + rhof*cos(phi5 + alphaf);
        y5 - y3 + rhof*sin(phi5 + alphaf);
        x6 - x4 + rhoe*cos(phi6 + alphae);
        y6 - y4 + rhoe*sin(phi6 + alphae);
        ];

PhiD=[phi5 + t;
      phi6 - t;
      x4 - l2*sin(phi4);
      y4 + l2*cos(phi4) - (l1 + l2)/2;
      ];

Phi=[PhiK;PhiD];

%Construção da jacobiana, nu e gamma

[Phiq, nu, Gamma] = construtor (q, qp, t, Phi);

%Os resultados tem que ser passados para as respectivas funções e nos nomes
%das variáveis convertidos para x(i).
%%
clear
%Entrada da estimativa da posição inicial e constantes do modelo
l1 = 5;
l2 = 5;
x0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;1;0;1;0];
t0 = 0;

%Início da simulação
x = x0;
t = t0;
i=2;

while t<2*pi
    %Análise de posição
    %Definir as funções funcaoseqnln e jacobiana
    A=jacobiana(x,t);
    if det(A)==0
    warning('jacobiana singular');
    end
    x=newtonraphson(x,t);

    %Análise de velocidade
    b=nu(x,t);
    A=jacobiana(x,t);
    xp = eliminacaogauss(A,b);

    %Análise de aceleração
    b = Gamma(x,xp,t);
    xpp = eliminacaogauss(A,b);

    %Atualização do tempo e estimativa da nova posição
    dt=0.0001;
    x = x + xp*dt + xpp*dt^2/2;
    t = t + dt;
    X(:,i)=[x;t];
    i=i+1;
end

%cames F
subplot(1,2,1)
polarplot(X(22,2:end),X(21,2:end))
title('Came F')

%cames E
subplot(1,2,2)
polarplot(X(20,2:end),X(19,2:end))
title('Came E')