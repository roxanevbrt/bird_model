classdef Stability < DefBalance
    %Code to calculate the balance once the no-load lengths have been selected

    properties
        dlePA=0;
        dleAr=0;
        dleAv=0;
        dleKH1=0;
        dleKH2=0;
    end

    methods

        function obj = Stability(dlePA,dleAr,dleAv,dleKH1,dleKH2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==5
                obj.dleAr=dleAr;
                obj.dleAv=dleAv;
                obj.dlePA=dlePA;
                obj.dleKH1=dleKH1;
                obj.dleKH2=dleKH2;
            end
        end

        function obj=defineOiseau(obj,L_f,L_tbt,L_tmt,rK1,rA,rP,L_h,dav,dar, Lg, Tg,m,rK2,Lfib)
               obj.L_f=L_f;
               obj.L_tmt=L_tmt;
               obj.L_tbt=L_tbt;
               obj.rK1=rK1;
               obj.rK2=rK2;
               obj.rA=rA;
               obj.rP=-rP;
               obj.doigt_lx2 = dav;
                obj.doigt_lx1 = -dar;
               obj.corps_joint_h=L_h;
               obj.L_g=Lg;
               obj.T_g=Tg*pi/180;
               obj.masse=m;
               obj.L_fibula=-Lfib;
        end



        function B=dU_modele_elastique(obj,q)
            %balance elastic potential term
            [l3,~] = eval_l4_12(obj,q);
            [~,dl4] = derivee_dlH_dqH(obj,q);
            dl4q3 = differencielle(obj, q,[0;0;1;0]);
            dl4q4 = differencielle(obj, q,[0;0;0;1]);
            if q(3)<pi
                theta=fminbnd(@(theta) lenght(obj,q,theta),q(3),pi);
                l42=lenght(obj,q,theta);
            else
                theta=fminbnd(@(theta) lenght(obj,q,theta),pi,q(3));
                l42=lenght(obj,q,theta);
            end
            lqPA=obj.rP*q(1)-obj.rA*q(2);
            lqAr=-obj.rA*q(2);
            lqAv=obj.rA/2*q(2);
            lqKH1=obj.rK1*q(3)+l3;
            lqKH2=l42;

            B=[obj.rP*obj.KPA*(lqPA-obj.dlePA);
                -(obj.KPA*(lqPA-obj.dlePA)+obj.KAr*(lqAr-obj.dleAr))*obj.rA+obj.rA/2*obj.KAv*(lqAv-obj.dleAv);
                obj.rK1*obj.KHK1*(lqKH1-obj.dleKH1)+dl4q3*obj.KHK2*(lqKH2-obj.dleKH2);
                +dl4*obj.KHK1*(lqKH1-obj.dleKH1)+dl4q4*obj.KHK2*(lqKH2-obj.dleKH2)];
        end
               
        function d=equilibre_elastique(obj,q)
            [dyP,dyA,dyK,dyH] = derivee_dy(obj,q); %gravity term
            BF=dU_modele_elastique(obj,q);
            A=obj.masse*obj.g*[dyP;dyA;dyK;dyH];
            d=A+BF; %balance equation

        end

        function qee= positionequilibre_elastique(obj,q)
            fsOpts = optimoptions('fsolve', 'MaxIter', 10000,'MaxFunEvals', 10000,'TolFun', 1e-15, 'TolX', 1e-6, 'Algorithm','levenberg-marquardt');
            qee=fsolve(@obj.equilibre_elastique,q,fsOpts);
        end

        function K = matrice_raideur_e(obj,q)
                %make the stability matrix
            Kg = matrice_Kg(obj,q);
            Ke=matrice_5cableKe_e(obj,q);
            K=Kg+Ke;
        end

        function K_e = matrice_5cableKe_e(obj,q)
            %Stability Matrix
            K_e1=zeros(4,4);
            K_e2r=zeros(4,4);
            K_e2v=zeros(4,4);
            K_e3=zeros(4,4);
            K_e4=zeros(4,4);

            K_e1(1,1)=obj.rP^2;
            K_e1(2,2)=obj.rA^2;
            K_e1(1,2)=-obj.rP*obj.rA;
            K_e1(2,1)=K_e1(1,2);

            K_e2v(2,2)=obj.rA^2/4;
            K_e2r(2,2)=obj.rA^2;

            [l4,dl4] = derivee_dlH_dqH(obj,q);
            K_e3(3,3)=obj.rK1^2;
            K_e3(4,4)=dl4^2;
            K_e3(3,4)=obj.rK1*dl4;
            K_e3(4,3)=K_e3(3,4);

            dl4q3 = differencielle(obj, q,[0;0;1;0]);
            dl4q4 = differencielle(obj, q,[0;0;0;1]);
            K_e4(3,3)=dl4q3^2;
            K_e4(4,4)=dl4q4^2;
            K_e4(3,4)=dl4q3*dl4q4;
            K_e4(4,3)=K_e4(3,4);

            K_e=obj.KPA*K_e1+obj.KAr*K_e2r+ obj.KAv*K_e2v+obj.KHK1*K_e3+ obj.KHK2*K_e4;

            eps=0.01;
            [dl44m, ~, dl444m]= eval_diff(obj,q-eps*[0;0;0;1]);
            [dl44p, ~, dl444p]=eval_diff(obj,q+eps*[0;0;0;1]);
            [~, dl443m]=eval_diff(obj,q-eps*[0;0;1;0]);
            [~, dl443p]=eval_diff(obj,q+eps*[0;0;1;0]);


            dl44_2=(dl44p-dl44m)/(2*eps);
            dl443_2=(dl443p-dl443m)/(2*eps);
            dl444_2=(dl444p-dl444m)/(2*eps);

            if q(3)<pi
                theta=fminbnd(@(theta) obj.lenght(q,theta),q(3),pi);
                l4_2=obj.lenght(q,theta);
            else
                theta=fminbnd(@(theta) obj.lenght(q,theta),pi,q(3));
                l4_2=obj.lenght(q,theta);
            end
            dl=[obj.rP*q(1)-obj.rA*q(2);-obj.rA*q(2);obj.rK1*q(3)+l4;l4_2];
            K_e(3,3)=K_e(3,3)+ obj.KHK2*dl443_2*(dl(4)-obj.dleKH2);
            K_e(4,4)=K_e(4,4)+ obj.KHK1*dl44_2*(dl(3)-obj.dleKH1)+obj.KHK2*dl444_2*(dl(4)-obj.dleKH2);
        end

        function Kg = matrice_Kg(obj,q)
            [d2y1, d2y2, d2y3, d2y4] = derivee_seconde(obj,q);
            % on construit la partie liee a la gravite
            Kg=[d2y1 d2y2 d2y3 d2y4 ; d2y2 d2y2 d2y3 d2y4 ; d2y3 d2y3 d2y3 d2y4 ; d2y4 d2y4 d2y4 d2y4];
            Kg =obj.masse*obj.g*Kg;
        end


        % Check Function
        function verification_stabilite_e(obj,q)
            K=matrice_raideur_e(obj,q);
            raid=eig(K);
            if min(raid)<0
                disp('The configuration is not stable with respect to a perturbation')
            else
                disp('The configuration is stable with respect to a perturbation')
            end
        end

        function a= etat_stabilite_e(obj,q)
            %function to test stability
            K=matrice_raideur_e(obj,q);
            raid=eig(K);
            if min(raid)<0
                a=0;
                %disp('The configuration is not stable with respect to a perturbation')
            else
                a=1;
                %disp('The configuration is stable with respect to a perturbation')
            end
        end


       function [tg, lg]=P_C(obj,q,CoM)
            CoMx=CoM(1);
            CoMy=CoM(2);
            P_hanche=obj.hip_position(q);
            Px=P_hanche(1);
            Py=P_hanche(2);

            syms TTg LLg
            eqn = [Px-CoMx-LLg*sin(q(1)+q(2)+q(3)+q(4)+TTg)==0,Py-CoMy+LLg*cos(q(1)+q(2)+q(3)+q(4)+TTg)==0,LLg>0];
            [NeTg, NeLg]=solve(eqn,[TTg LLg]);

            tg=double(NeTg);
            lg=double(NeLg);
        end

    end
end