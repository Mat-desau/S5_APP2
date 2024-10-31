clc
clear all
close all

%% GLISSADE

% Calcul du polynome

E = 12.5;

xn = [0 8 15 20 25];
yn = [30 19 20 16 E];

C =  [ 1 0 0 0 0;
       1 8 64 512 4096;
       1 15 225 3375 50625;
       1 20 400 8000 160000;
       1 25 625 15625 390625]

S = [30;
     19;
     20;
     16;
     E]

A = pinv(C)*S;

H = @(x) A(5)*x.^4 + A(4)*x.^3 + A(3)*x.^2 + A(2)*x + A(1);

% Affichage courbe approximée 

fplot(H, [0 25], 'LineWidth', 2, 'Color', 'red')
hold on
scatter(xn, yn, "Filled", "red")
xlim([0 30])
ylim([10 32])
title("Courbe approximée par un polynôme")
xlabel("Distance (m)")
ylabel("Hauteur (m)")


% Dérivée du polynome

Hd = @(x) 4*A(5)*x.^3 + 3*A(4)*x.^2 + 2*A(3)*x + A(2);

% Trouver extremums de la fonction
% Racines de la dérivée

Ext = roots([4*A(5) 3*A(4) 2*A(3) A(2)])

% Trouver les y associées aux extremums

y1 = A(5)*Ext(1).^4 + A(4)*Ext(1).^3 + A(3)*Ext(1).^2 + A(2)*Ext(1) + A(1)
y2 = A(5)*Ext(2).^4 + A(4)*Ext(2).^3 + A(3)*Ext(2).^2 + A(2)*Ext(2) + A(1)
y3 = A(5)*Ext(3).^4 + A(4)*Ext(3).^3 + A(3)*Ext(3).^2 + A(2)*Ext(3) + A(1)
yE = [y1 y2 y3];

% Affichage des extremums

figure
fplot(H, [0 25], 'LineWidth', 2, 'Color', 'black')
hold on
scatter(Ext, yE, "Filled", "red")
xlim([0 30])
ylim([10 32])
title("Extremums sur la glissade")
xlabel("Distance (m)")
ylabel("Hauteur (m)")


% Trouver les coefficients de frottements max et min

yi = 30;
g = 9.81;
vmax45 = 45/3.6;
vmax25 = 25/3.6;
vmin20 = 20/3.6;
vmin15 = 15/3.6;
vmin10 = 10/3.6;
xmax = 25;

% Au premier max de vitesse
umin1 = ((g*(yi-y3))-(1/2*vmax45.^2))/(Ext(3)*g)

% Au minium de vitesse
umax = ((g*(yi-y2))-(1/2*vmin10.^2))/(Ext(2)*g)

% Au point E max de vitesse 25kmh
uminE = ((g*(yi-y1))-(1/2*vmax25.^2))/(Ext(1)*g)

% Au point E min de vitesse 20kmh
umaxE = ((g*(yi-y1))-(1/2*vmin20.^2))/(Ext(1)*g)

% uc choisi
uc = 0.62

% uc avec ERMS
ucmax = uc + 0.018
ucmin = uc - 0.018

% Graphique vitesse

V = @(x)sqrt((yi-uc*x-(A(5)*x.^4 + A(4)*x.^3 + A(3)*x.^2 + A(2)*x + A(1)))*2*g);

% Calcul de vitesse avec erreur ERMS

Vfin = sqrt((yi-uc*xmax-(A(5)*xmax.^4 + A(4)*xmax.^3 + A(3)*xmax.^2 + A(2)*xmax + A(1)))*2*g)
Vfin_min = sqrt((yi-ucmax*xmax-(A(5)*xmax.^4 + A(4)*xmax.^3 + A(3)*xmax.^2 + A(2)*xmax + A(1)))*2*g)
Vfin_max = sqrt((yi-ucmin*xmax-(A(5)*xmax.^4 + A(4)*xmax.^3 + A(3)*xmax.^2 + A(2)*xmax + A(1)))*2*g)

Vfin_kmh = Vfin * 3.6
Vfin_kmh_min = Vfin_min * 3.6
Vfin_kmh_max = Vfin_max * 3.6

figure
fplot(V, [0 25], 'LineWidth', 2, 'Color', 'blue')
xlim([0 30])
ylim([0 16])
title("Vitesse de l'usager")
xlabel("Distance (m)")
ylabel("Vitesse (m/s)")

%% POLYNOME DE LISSAGE DE LA VALVE

uc = 0.62;
Pourc = [0 10 20 30 40 50 60 70 80 90 100];
uct = [0.87 0.78 0.71 0.61 0.62 00.51 0.51 0.49 0.46 0.48 0.46];

C = [1 0 0;
     1 10 100;
     1 20 400;
     1 30 900;
     1 40 1600;
     1 50 2500;
     1 60 3600;
     1 70 4900;
     1 80 6400;
     1 90 8100;
     1 100 10000]

S = [0.87;
     0.78;
     0.71;
     0.61;
     0.62;
     0.51;
     0.51;
     0.49;
     0.46;
     0.48;
     0.46]

format long

A = pinv(C)*S


F = @(x) A(3)*x.^2 + A(2)*x + A(1);

fplot(F, 'LineWidth', 2, 'Color', 'red')
hold on
scatter(Pourc, uct, "Filled", "red")
xlim([0 100])
ylim([0 1])
title("Courbe de lissage des données de frottement")
xlabel("Ouverture de la valve (%)")
ylabel("Frottement μ_c")

% Calcul du ERMS

N = 11;

g = A(3)*Pourc.^2+A(2)*Pourc+A(1)

ERMS = sqrt(1/N*sum((g-uc).^2))

RacF = roots([A(3) A(2) A(1)-uc])

% 33.19% d'ouverture de la valve

%% MINUTERIE