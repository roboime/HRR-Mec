% fun��o x = 2*l1 y = - 2*l2
%resultado conhecido = 2 c�rculos
%escolha dos valores de forma a evitar singularidades
% w4 = 1
%%Simulador de mecanismos planos
clear
%%
%Constru��o do modelo baseado em vari�veis simb�licas
%Defini��o das vari�veis
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4 rho1 alpha1 rho2 alpha2 real %posi��o
syms x1p y1p phi1p x2p y2p phi2p x3p y3p phi3p x4p y4p phi4p rho1p alpha1p rho2p alpha2p real %velocidade
syms l1 l2 %constantes do modelo
syms t %tempo

%Constru��o do vetor de coordenadas generalizadas e a sua derivada
q = [x1;y1;phi1;x2;y2;phi2;x3;y3;phi3;x4;y4;phi4;rho1;alpha1;rho2;alpha2];

qp = [x1p;y1p;phi1p;x2p;y2p;phi2p;x3p;y3p;phi3p;x4p;y4p;phi4p;rho1p;alpha1p;rho2p;alpha2p];

%Constru��o do vetor de restri��es cinem�ticas e diretoras
PhiK = [x4;
        y4;
        y1;
        x2;
        phi1;
        phi2;
        phi3;
        y3 - y2;
        x3 - x1;
        x4 - x1 + rho1*cos(phi4 + alpha1) + l1*cos(phi1);
        y4 - y1 + rho1*sin(phi4 + alpha1) + l1*sin(phi1);
        x4 - x2 + rho2*cos(phi4 + alpha2) + l2*sin(phi2);
        y4 - y2 + rho2*sin(phi4 + alpha2) - l2*cos(phi2)
        ];

PhiD=[x3 - 2*l1;
      y3 + 2*l2;
      phi4 - t
      ];

Phi=[PhiK;PhiD];

%Constru��o da jacobiana, nu e gamma

[Phiq, nu, Gamma] = construtor (q, qp, t, Phi);

%Os resultados tem que ser passados para as respectivas fun��es e nos nomes
%das vari�veis convertidos para x(i).
%%
clear
%Entrada da estimativa da posi��o inicial e constantes do modelo
l1 = 5;
l2 = 3;
x0 = [0;0;0;0.4;0;1;0;0;0;0;0.7;0;1;0;1;0];
t0 = 0;

%In�cio da simula��o
x = x0;
t = t0;
i=2;

while t<2*pi
    %An�lise de posi��o
    %Definir as fun��es funcaoseqnln e jacobiana
    A=jacobiana(x,t);
    if det(A)==0
    warning('jacobiana singular');
    end
    x=newtonraphson(x,t);

    %An�lise de velocidade
    b=nu(x,t);
    A=jacobiana(x,t);
    xp = eliminacaogauss(A,b);

    %An�lise de acelera��o
    b = Gamma(x,xp,t);
    xpp = eliminacaogauss(A,b);

    %Atualiza��o do tempo e estimativa da nova posi��o
    dt=0.0001;
    x = x + xp*dt + xpp*dt^2/2;
    t = t + dt;
    X(:,i)=[x;t];
    i=i+1;
end

polarplot(X(14,2:end),X(13,2:end))
hold on
polarplot(X(16,2:end),X(15,2:end)) 
legend('Came 1','Came 2')
hold off