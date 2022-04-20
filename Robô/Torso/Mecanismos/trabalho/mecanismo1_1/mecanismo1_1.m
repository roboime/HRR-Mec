%fun��o x2 = l2 + ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*cos(t) e y2 = - l1 + ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*sin(t)
%fun��o escolhida para n�o haver singularidades
%%Simulador de mecanismos planos
clear
%%
%Constru��o do modelo baseado em vari�veis simb�licas
%Defini��o das vari�veis
syms x1 y1 phi1 x2 y2 phi2 real %posi��o
syms x1p y1p phi1p x2p y2p phi2p real %velocidade
syms l l1 l2 %constantes do modelo
syms t %tempo

%Constru��o do vetor de coordenadas generalizadas e a sua derivada
q = [x1;y1;phi1;x2;y2;phi2];

qp = [x1p;y1p;phi1p;x2p;y2p;phi2p];

%Constru��o do vetor de restri��es cinem�ticas e diretoras
PhiK = [x1;
        y1;
        x1 + l1*cos(phi1) - x2 + l2*cos(phi2)
        y1 + l1*sin(phi1) - y2 + l2*sin(phi2)];

PhiD=[x2 - l2 - ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*cos(t);
      y2 + l1 - ((l1 + l2 - sqrt(l1^2 + l2^2))/2)*sin(t)];

Phi=[PhiK;PhiD];

%Constru��o da jacobiana, nu e gamma

[Phiq, nu, Gamma] = construtor (q, qp, t, Phi);

%Os resultados tem que ser passados para as respectivas fun��es e nos nomes
%das vari�veis convertidos para x(i).
%%
clear
%Entrada da estimativa da posi��o inicial e constantes do modelo
l1 = 5;
l2 = 5;
l = 10;
x0 = [0;0;0.1;0;0;0];
t0 = 0;

%In�cio da simula��o
x = x0;
t = t0;
V=[];
Y=[];
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
    
    %velocidade do cursor perpendicular � barra
    vc = -xp(4)*sin(x(6))+xp(5)*cos(x(6))-l*xp(6);
    V(:,i) = vc;

    %Atualiza��o do tempo e estimativa da nova posi��o
    dt=0.0001;
    x = x + xp*dt + xpp*dt^2/2;
    t = t + dt;
    X(:,i)=[x;t];
    i=i+1;
end

%�rea da fun��o
fprintf('A �rea da curva � %.4f\n',polyarea(X(4,2:end),X(5,2:end)))

%lrtheta
fprintf('l*r*theta = %.4f\n',l2*trapz(X(7,2:end),V(2:end)))
