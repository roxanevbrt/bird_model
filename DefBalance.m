
classdef DefBalance < DefCoM
% define free lenght and force ratio

    properties
        KPA=70;
        KAr=70;
        KAv=70;
        KHK1=15;
        KHK2=15;

    end

    methods
        function obj = DefBalance(KPA,KAr,KAv, KHK1,KHK2)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin ==5
                obj.KPA = KPA;
                obj.KAr = KAr;
                obj.KAv = KAv;
                obj.KHK1 = KHK1;
                obj.KHK2 = KHK2;
            end
        end

               function [F1, F2, F4, F5]=rapportforce(obj,q,F3)
           %function to find the needed force to reach balance
            [~,dl4] = derivee_dlH_dqH(obj,q);
            dl4q3 = differencielle(obj, q,[0;0;1;0]);
            dl4q4 = differencielle(obj, q,[0;0;0;1]);
            [dy1,dy2,dy3,dy4]=obj.derivee_dy(q);

            F4=(dy1*obj.g*obj.masse)/abs(obj.rP);
            alpha=(dy2*obj.g*obj.masse)/obj.rA-F4;   %to define f3  +F5/2;
            F5=(F3-alpha)*2;
            syms F1b F2b
            eqn = [dy3*obj.g*obj.masse+obj.rK1*F1b+dl4q3*F2b==0,dy4*obj.g*obj.masse+dl4*F1b+dl4q4*F2b==0,F1b>0,F2b>0];
            [F3s, F4s]=solve(eqn,[F1b F2b]);
            F1=double(F3s);
            F2=double(F4s);
               end

         function [dl1,dl2,dl3,dl4,dl5] = deflongueuravide(obj,q,F_ini)
            [l4,~] = eval_l4_12(obj,q);
            theta=fminbnd(@(theta) lenght(obj,q,theta),q(3),pi);
            l4_2=lenght(obj,q,theta);
            dl_ini1=obj.rK1*q(3)+l4;
            dl_ini2=l4_2;
            dl_ini3=-obj.rA*q(2);
            dl_ini4=obj.rP*q(1)-obj.rA*q(2);

            dl1=dl_ini1-F_ini(1)/obj.KHK1;
            dl2=dl_ini2-F_ini(2)/obj.KHK2;
            dl3=dl_ini3-F_ini(3)/obj.KAr;
            dl4=dl_ini4-F_ini(4)/obj.KPA;
            if F_ini(5)==0
                dl5=0;
            else
                dl_ini5=obj.rA/2*q(2);
                dl5=dl_ini5-F_ini(5)/obj.KAv;
            end
            

        end
        
        function verification_force(obj,Forces)
            if min(Forces)>=0
                disp('The forces in the cablee are positive, the value of F/mg is :')
                disp(Forces/(obj.masse*obj.g))
            else
                msg = "The force is the cable is negativ, the equilibrium is impossible";
                msgbox(msg);
                error(msg)
            end
        end

        function a = etat_force(obj,Forces)
            if min(Forces)>=0
                a=1;
                %disp(Forces/(obj.masse*obj.g))
            else
                a=0;
            end
        end

    end
end