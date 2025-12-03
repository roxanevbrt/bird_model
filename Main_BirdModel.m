%%Programme to test the feasibility of balance and mobility of the centre of mass of different species

clear

%To change the model, uncomment the associated bird and choose the correct model for the database (l.39).
%Zebra finch values 
% q=[0.89; -1.4158; 1.65; -2.1489]; %Desired position
% Stif_max=16; %Stiffness choice for K0
% limite=10; %Variation of x and y CoM
% Variation_x=-limite:1:limite;
% Variation_y=-limite:1:limite;

%Malard values 
% q=[0.13; -0.66; 1.92; -2.1489]; 
% Stif_max=50;
% limite=50;
% Variation_x=-65:5:limite;
% Variation_y=-limite:5:limite;
                
%Egret values
% q=[0.087; -0.29; 1.98; -2.1489]; %Egret position
% Stif_max=100;
% limite=50;
% Variation_x=-70:10:limite;
% Variation_y=-limite:10:limite;

%Military macaw values
q=[1.12; -1.68; 1.83; -1.7]; 
Stif_max=800;
limite=40;
Variation_x=-55:5:limite;
Variation_y=-limite:5:limite;


qe=q;
MonOiseau=Stability;
% 1-Zebra finch, 2-Egret, 3-Malard, 4-Military macaw
[L_f,L_tbt,L_tmt,rK1,rA,rP,L_h,name, masse, rK2, Lfib, dav,dar, Lg, Tg] = modelchoice(4); %change for good model choice
MonOiseau=MonOiseau.defineOiseau(L_f,L_tbt,L_tmt,rK1,rA,rP,L_h,dav,dar, Lg, Tg,masse,rK2,Lfib);


% CoM Position
CoM=MonOiseau.CoM(qe);

figure() 
MonOiseau.dessin_patte(qe);
title(['Desired equilibrium position of '  name])
grid off

F3=1; %Choice for F3
[F1, F2, F4, F5]=MonOiseau.rapportforce(qe,F3); %force ratio to obtain equilibrium position
Forces=[F1,F2,F3,F4,F5];
Forces_ration=Forces/(MonOiseau.masse*MonOiseau.g);
     
%Stiffness choice (16 Zebra finch; 50 Malard ; 100 Egret ; 800 Military
%macaw)
MonOiseau.KPA=Stif_max; %stiffness for C4
MonOiseau.KAr=Stif_max/2; %stiffness for C3
MonOiseau.KAv=Stif_max/4; %stiffness for C5 =0 for 4cables model
MonOiseau.KHK1=Stif_max/8; %stiffness for C1
MonOiseau.KHK2=Stif_max/8; %stiffness for C2

% We calculate free length using the forces we found.
[dle1,dle2,dle3,dle4,dle5] = MonOiseau.deflongueuravide(qe,Forces);
%Affectation of free length
MonOiseau.dlePA=dle4;
MonOiseau.dleAr=dle3;
MonOiseau.dleAv=dle5;
MonOiseau.dleKH1=dle1;
MonOiseau.dleKH2=dle2;

%Check if it's the desired position 
qee=MonOiseau.positionequilibre_elastique(qe);
figure()
MonOiseau.dessin_patte(qee);
title(['Equilibrium position of '  name])

%Check forces positivity 
MonOiseau.verification_force(Forces);
%On verifie la stabilité mecanique
MonOiseau.verification_CoM(qee);
%On verifie la stabilite automatique
MonOiseau.verification_stabilite_e(qee);


%% Variation of CoM
%Stability check
CoMte=zeros(length(Variation_x),length(Variation_y));    %Mecanique
ee=zeros(length(Variation_x),length(Variation_y));       %Automatique
aa=zeros(length(Variation_x),length(Variation_y));  
%Position CoM 
CoMix=length(Variation_x);
CoMiy=length(Variation_y);


figure() %figure at previous equilibrium before changing the CoM
MonOiseau.dessin_patte_point(qe)
hold on
ax_limits=axis;
PCOMtotxx=zeros(length(Variation_x),length(Variation_y));

for i=1:1:length(Variation_x)
    for j=1:1:length(Variation_y)
        CoMi=CoM+[Variation_x(i), Variation_y(j)];
        CoMix(i)=CoM(1)+Variation_x(i);
        CoMiy(j)=CoM(2)+Variation_y(j);
        [Tgi,Lgi]=MonOiseau.P_C(qe,CoMi);
        MonOiseau.T_g=Tgi;
        MonOiseau.L_g=Lgi;
        qee=MonOiseau.positionequilibre_elastique(qe);
        CoMte(i,j)=MonOiseau.etat_CoM(qee);
        %check stability
        ee(i,j)=MonOiseau.etat_stabilite_e(qee);
        PCOM=MonOiseau.CoM(qee);
        PCOMtotxx(i,j)=PCOM(1);
    end
end

v=[0,0];
v2=[MonOiseau.doigt_lx2,MonOiseau.doigt_lx2];
v3=[MonOiseau.doigt_lx1,MonOiseau.doigt_lx1];
contour(CoMix, CoMiy,PCOMtotxx',v,'k','LineWidth',1)
hold on
contour(CoMix, CoMiy,PCOMtotxx',v2,'k','LineWidth',1)
hold on
contour(CoMix, CoMiy,PCOMtotxx',v3,'k','LineWidth',1)
hold on
plot([CoMix(1) CoMix(1) CoMix(length(Variation_x)) CoMix(length(Variation_x)) CoMix(1)],[CoMiy(1) CoMiy(length(Variation_y)) CoMiy(length(Variation_y)) CoMiy(1) CoMiy(1) ],'k','LineWidth',1);
ax = gca;
ax.XLim=ax_limits(1:2);

% Extraction des coordonnées des contours
[C1, h1] = contour(CoMix, CoMiy, PCOMtotxx', v);
[C2, h2] = contour(CoMix, CoMiy, PCOMtotxx', v2);
[C3, h3] = contour(CoMix, CoMiy, PCOMtotxx', v3);

% Extraction des données des contours
x1 = C1(1, 2:end); y1 = C1(2, 2:end);
x2 = C2(1, 2:end); y2 = C2(2, 2:end);
x3 = C3(1, 2:end); y3 = C3(2, 2:end);
x_front_fall=[ CoMix(1) CoMix(length(Variation_x)) CoMix(length(Variation_x)) fliplr(x2)];
y_front_fall = [CoMiy(length(Variation_y)) CoMiy(length(Variation_y)) CoMiy(1) fliplr(y2) ];
x_rear_fall=[ CoMix(1) fliplr(x3) CoMix(1) CoMix(1) ];
y_rear_fall = [CoMiy(1) fliplr(y3) CoMiy(length(Variation_y)) CoMiy(1) ];
x_good=[x1, fliplr(x2)];
y_good= [y1, fliplr(y2)];

hold on;
fill(x_good, y_good, [0, 0.8, 0], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % green zone
fill([x3, fliplr(x1)], [y3, fliplr(y1)],[0.6 0 .6] , 'FaceAlpha', 0.7, 'EdgeColor', 'none'); %purple zone
fill(x_front_fall, y_front_fall, [0 .2 .7], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
fill(x_rear_fall, y_rear_fall, [0.6, 0.8, 1], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
grid on;
hold off;
title(['Balance position with K1 =' num2str(MonOiseau.KHK1) ', K2=' num2str(MonOiseau.KHK2) ', K3=' num2str(MonOiseau.KAr) ', K4=' num2str(MonOiseau.KPA) , ', K5=' num2str(MonOiseau.KAv) ])

