classdef CaracteristiqueOiseau
    %define bird model

    properties
        L_f=237.24;
        L_tbt=451.72;
        L_tmt=328.74;
        L_fibula=-2*1.8*50;
        rK1= 50 ;
        rK2 = 50;
        rA= 40 ;
        rP=-20;
        corps_joint_h=97;
        corps_joint_l=0;

    end

    methods

        function obj = CaracteristiqueOiseau(L_f,L_tbt,L_tmt,rK1,rK2,rA,rP,corps_joint_h,corps_joint_l)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 9
                obj.L_f = L_f;
                obj.L_tbt = L_tbt;
                obj.L_tmt = L_tmt;
                obj.rK1 = rK1;
                obj.rK2 = rK2;
                obj.rA = rA;
                obj.rP = rP;
                obj.corps_joint_h = corps_joint_h;
                obj.corps_joint_l = corps_joint_l;
            end
        end

        function f=lenght(obj,q,theta)
            %Find position of ligament loop
            % The points are expressed in a reference frame linked to the
            % tibia with the knee as the origin
            P10=[0; obj.L_fibula];
            Dist=sqrt(obj.L_f^2+obj.corps_joint_h^2-2*obj.L_f*obj.corps_joint_h*cos(q(4)));
            alpha= q(3)+atan2(-obj.corps_joint_h*sin(q(4)),obj.L_f-obj.corps_joint_h*cos(q(4)));
            P7=Dist*[-sin(alpha);cos(alpha)];
            P11=obj.rK2*[-sin(theta);cos(theta)];
            f=norm(P7-P11)+norm(P10-P11);
        end

        function [l4,dl4] = derivee_dlH_dqH(obj,q)
            %dlH(qH)/dqH
            V=-obj.corps_joint_h*sin(q(4))+obj.corps_joint_l*cos(q(4));
            U=obj.L_f-obj.corps_joint_h*cos(q(4))-obj.corps_joint_l*sin(q(4));
            la=sqrt(obj.L_f^2+obj.corps_joint_h^2+obj.corps_joint_l^2-2*obj.L_f*obj.corps_joint_h*cos(q(4))-2*obj.L_f*obj.corps_joint_l*sin(q(4)));
            a41=atan2(U,V);
            dU=obj.corps_joint_h*sin(q(4))-obj.corps_joint_l*cos(q(4));
            dV=-obj.corps_joint_h*cos(q(4))-obj.corps_joint_l*sin(q(4));
            dla=(2*obj.L_f*obj.corps_joint_h*sin(q(4))-2*obj.L_f*obj.corps_joint_l*cos(q(4)))/la/2;
            da41=-(dV*U-dU*V)/(U*U+V*V);
            
            a42=-acos(obj.rK1/la);
            a4=a41+a42;
            l4=sqrt(la^2-obj.rK1^2)+obj.rK1*a4;
            
            da42=-obj.rK1*dla/(la*sqrt(la*la-obj.rK1*obj.rK1));
            da4=da41+da42;
            dl4=la/(sqrt(la^2-obj.rK1^2))*dla+obj.rK1*da4;
        end

        function dl4q = differencielle(obj, q,e)
            eps=0.001;
            e=eps*e;
            if q(3)<pi
                thetam=fminbnd(@(theta) obj.lenght(q-e,theta),q(3),pi);
                lm=obj.lenght(q-e,thetam);
                thetap=fminbnd(@(theta) obj.lenght(q+e,theta),q(3),pi);
                lp=obj.lenght(q+e,thetap);
            else 
                thetam=fminbnd(@(theta) obj.lenght(q-e,theta),pi,q(3));
                lm=obj.lenght(q-e,thetam);
                thetap=fminbnd(@(theta) obj.lenght(q+e,theta),pi,q(3));
                lp=obj.lenght(q+e,thetap);
            end
            dl4q=(lp-lm)/(2*eps);
        end

        function [dl4,dl4q3, dl4q4]= eval_diff(obj,q)
            [~,dl4] = derivee_dlH_dqH(obj,q);
            dl4q3 = differencielle(obj, q,[0;0;1;0]);
            dl4q4 = differencielle(obj, q,[0;0;0;1]);
        end

        function [l4,l4_2]  = eval_l4_12(obj,q)
            q4=q(4);
            V=-obj.corps_joint_h*sin(q4)+obj.corps_joint_l*cos(q4);
            U=obj.L_f-obj.corps_joint_h*cos(q4)-obj.corps_joint_l*sin(q4);
            la=sqrt(obj.L_f^2+obj.corps_joint_h^2+obj.corps_joint_l^2-2*obj.L_f*obj.corps_joint_h*cos(q4)-2*obj.L_f*obj.corps_joint_l*sin(q4));
            a41=atan2(U,V);

            %first cable
            a42=-acos(obj.rK1/la);
            a4=a41+a42;
            l4=sqrt(la^2-obj.rK1^2)+obj.rK1*a4;

            %second cable 
            a42_2=-acos(obj.rK2/la);
            a4_2=a41+a42_2;
            l4_2=sqrt(la^2-obj.rK2^2)+obj.rK2*a4_2;
        end

        function [dyP,dyA,dyK,dyH] = derivee_dy(obj,q)
            dyP=-obj.L_tmt*sin(q(1))-obj.L_tbt*sin(q(1)+q(2))-obj.L_f*sin(q(1)+q(2)+q(3))-obj.L_g*sin(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            dyA=-obj.L_tbt*sin(q(1)+q(2))-obj.L_f*sin(q(1)+q(2)+q(3))-obj.L_g*sin(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            dyK=-obj.L_f*sin(q(1)+q(2)+q(3))-obj.L_g*sin(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            dyH=-obj.L_g*sin(q(1)+q(2)+q(3)+q(4)+obj.T_g);
        end

        function [d2y1, d2y2, d2y3, d2y4] = derivee_seconde(obj,q)
            d2y1=-obj.L_tmt*cos(q(1))-obj.L_tbt*cos(q(1)+q(2))-obj.L_f*cos(q(1)+q(2)+q(3))+...
                -obj.L_g*cos(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            d2y2=-obj.L_tbt*cos(q(1)+q(2))-obj.L_f*cos(q(1)+q(2)+q(3))+...
                -obj.L_g*cos(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            d2y3=-obj.L_f*cos(q(1)+q(2)+q(3))+...
                -obj.L_g*cos(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            d2y4= -obj.L_g*cos(q(1)+q(2)+q(3)+q(4)+obj.T_g);
        end


    
 
    
    end
end