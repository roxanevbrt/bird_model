classdef DefCoM < CaracteristiqueOiseau
    %make figure

    properties
        L_g=285;
        T_g=-45*pi/180;
        masse=0.134;
        g=9.81;
        doigt_lx1 =-50;
        doigt_lx2=100;
        doigt_ly =20;
    end

    methods
        function obj = DefCoM(L_g,T_g,masse,g,doigt_lx1, doigt_lx2,doigt_ly)
            if nargin == 7
                obj.L_g = L_g;
                obj.T_g = T_g;
                obj.masse = masse;
                obj.g = g;
                obj.doigt_lx1 = doigt_lx1;
                obj.doigt_lx2 = doigt_lx2;
                obj.doigt_ly = doigt_ly;
            end
        end

        function P_CoM=CoM(obj,q)
            %function to find the center of masse position
            x_CoM= -obj.L_tmt*sin(q(1))-obj.L_tbt*sin(q(1)+q(2))-obj.L_f*sin(q(1)+q(2)+q(3))-obj.L_g*sin(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            y_CoM= abs(obj.rP)+obj.L_tmt*cos(q(1))+obj.L_tbt*cos(q(1)+q(2))+obj.L_f*cos(q(1)+q(2)+q(3))+obj.L_g*cos(q(1)+q(2)+q(3)+q(4)+obj.T_g);
            P_CoM=[x_CoM y_CoM];
        end

        function P_hanche=hip_position(obj,q)
            x_hip= -obj.L_tmt*sin(q(1))-obj.L_tbt*sin(q(1)+q(2))-obj.L_f*sin(q(1)+q(2)+q(3));
            y_hip= abs(obj.rP)+obj.L_tmt*cos(q(1))+obj.L_tbt*cos(q(1)+q(2))+obj.L_f*cos(q(1)+q(2)+q(3));
            P_hanche=[x_hip y_hip];
        end

        function verification_CoM (obj, q)
            P_CoM= obj.CoM(q);
            [x_CoM]=P_CoM(1);
            if x_CoM<=obj.doigt_lx2  && x_CoM>=obj.doigt_lx1
                disp('The bird is in equilibrium on its fingers')
            elseif x_CoM>obj.doigt_lx2
                msg = "The bird tilts forward";
                f = msgbox(msg);
                error(msg)
            elseif x_CoM<obj.doigt_lx1
                msg = "The bird tilts backward";
                f = msgbox(msg);
                error(msg)
            end
        end

        function e=etat_CoM (obj, q)
            P_CoM= obj.CoM(q);
            [x_CoM]=P_CoM(1);
            if x_CoM<=obj.doigt_lx2  && x_CoM>=obj.doigt_lx1
                e=1;
                %disp('The bird is in equilibrium on its fingers')
            elseif x_CoM>obj.doigt_lx2
                e=0;
            elseif x_CoM<obj.doigt_lx1
                e=3;
            end
        end
        
        function dessin_patte(obj,q)
            %Foot draw
            P3=[0;abs(obj.rP)];
            plot([obj.doigt_lx1 obj.doigt_lx2 P3(1) obj.doigt_lx1], [0 0 P3(2) 0],'k','LineWidth',2)
            hold on

            %joint position
            P_cheville= P3 + obj.L_tmt*[-sin(q(1)); cos(q(1))];
            P_genou= P_cheville + obj.L_tbt*[-sin(q(1)+q(2)); cos(q(1)+q(2))];
            P_hanche= P_genou + obj.L_f*[-sin(q(1)+q(2)+q(3)); cos(q(1)+q(2)+q(3))];
            plot([P3(1) P_cheville(1) P_genou(1) P_hanche(1)],[P3(2) P_cheville(2) P_genou(2) P_hanche(2)],'k','LineWidth',2)
            %pulley
            plotcircle(P_cheville,obj.rA)
            plotcircle(P_genou,obj.rK1)
            plotcircle(P3,obj.rP)
            plotcirclepoint(P_cheville,obj.rA/2)
            grid

            obj.plotcorps(P_hanche,q)
            obj.plotCoM(q)

            % yellow cable
            if obj.rP>0
                plot([obj.rP obj.rP+obj.doigt_lx2/3],[obj.rP 2*obj.rP/3],'color',[255/255 164/255 36/255],'LineWidth',2)
                beta = 0;
            else
                plot([0 abs(obj.rP)+obj.doigt_lx2/3],[0 0],'color',[255/255 164/255 36/255],'LineWidth',2)
                beta=pi/2;
            end
            lo=sqrt((obj.L_tmt)^2-(obj.rP+obj.rA)^2);
            beta_2=-asin((obj.rA+obj.rP)/obj.L_tmt);
            R1=(obj.rA*P3+obj.rP*P_cheville)/(obj.rP+obj.rA)-obj.rP/(obj.rP+obj.rA)*[lo*sin(-q(1)+beta_2);lo*cos(-q(1)+beta_2)];
            R2=(obj.rA*P3+obj.rP*P_cheville)/(obj.rP+obj.rA)+obj.rA/(obj.rP+obj.rA)*[lo*sin(-q(1)+beta_2);lo*cos(-q(1)+beta_2)];
            plot([R1(1) R2(1)],[R1(2) R2(2)],'color',[255/255 164/255 36/255],'LineWidth',2)
            a=beta:(q(1)-beta_2-beta)/100:q(1)-beta_2;
            plot(P3(1)+obj.rP*cos(a),P3(2)+obj.rP*sin(a),'color',[255/255 164/255 36/255],'LineWidth',2); %jonction segment au niveau du pied
            beta=asin((obj.rK1)/(3/4*obj.L_tbt));
            R3=P_cheville-obj.rA*[cos(q(1)+q(2)-beta);sin(q(1)+q(2)-beta)];
            R4y=(P_cheville+3*P_genou)/4;
            plot([R3(1) R4y(1)],[R3(2) R4y(2)],'color',[255/255 164/255 36/255],'LineWidth',2)
            a=q(1)-beta_2:(q(2)-beta+beta_2)/100:q(1)+q(2)-beta;

            %green cable
            hv=0.9;
            beta=asin((obj.rA)/(hv*obj.L_tbt));
            R3=P_cheville-obj.rA*[cos(q(1)+q(2)-beta);sin(q(1)+q(2)-beta)];
            R4=((1-hv)*P_cheville+hv*P_genou);
            plot([R3(1) R4(1)],[R3(2) R4(2)],'color',[0 145/255 69/255],'LineWidth',2)
            hp=0.6;
            beta2=-asin((obj.rA)/(hp*obj.L_tmt));
            R3=P_cheville-obj.rA*[cos(q(1)-beta2);sin(q(1)-beta2)];
            R4=((1-hp)*P3+hp*P_cheville);
            plot([R3(1) R4(1)],[R3(2) R4(2)],'color',[0 145/255 69/255],'LineWidth',2)
            a=q(1)-beta2:(q(2)-beta+beta2)/100:q(1)+q(2)-beta;
            plot(P_cheville(1)- obj.rA*cos(a),P_cheville(2)- obj.rA*sin(a),':','color',[0 145/255 69/255],'LineWidth',2);
    
            %blue cable
            hv=0.9;
            beta=asin((-obj.rA/2)/(hv*obj.L_tbt));
            R3=P_cheville+obj.rA/2*[cos(q(1)+q(2)-beta);sin(q(1)+q(2)-beta)];
            R4=((1-hv)*P_cheville+hv*P_genou);
            plot([R3(1) R4(1)],[R3(2) R4(2)],'color',[65/255, 105/255, 225/255],'LineWidth',2)
            a=-pi()/2:pi()/180:q(1)+q(2)-beta;
            plot(P_cheville(1)+ obj.rA/2*cos(a),P_cheville(2)+ obj.rA/2*sin(a),'color',[65/255, 105/255, 225/255],'LineWidth',2);

            %red cable
            Dist=sqrt(obj.L_f^2+obj.corps_joint_h^2-2*obj.L_f*obj.corps_joint_h*cos(q(4)));
            alpha=q(1)+q(2)+q(3)+atan2(-obj.corps_joint_h*sin(q(4)),obj.L_f-obj.corps_joint_h*cos(q(4)));
            P7=P_genou+Dist*[-sin(alpha);cos(alpha)];
            a4=+atan2(-obj.corps_joint_h*sin(q(4)),obj.L_f-obj.corps_joint_h*cos(q(4)))-acos(obj.rK1/Dist);
            gamma2=pi/2+ q(1)+q(2)+q(3)+a4;
            a=q(1)+q(2):abs((gamma2-q(1)+q(2)))/100:gamma2;
            plot(P_genou(1)+ obj.rK1*cos(a),P_genou(2)+ obj.rK1*sin(a),'color',[1 0 46/255],'LineWidth',2);
            R8=P_genou+obj.rK1*[cos(gamma2);sin(gamma2)];
            plot([P7(1) R8(1)],[P7(2) R8(2)],'color',[1 0 46/255],'LineWidth',2)

            %on dessine le ligament rose
            P10=P_cheville + (obj.L_tbt+obj.L_fibula)*[-sin(q(1)+q(2)); cos(q(1)+q(2))];
            theta=fminbnd(@(theta) obj.lenght(q,theta),q(3),pi);
            P11=P_genou+obj.rK2*[-sin((q(1)+q(2)+theta));cos((q(1)+q(2)+theta))];
            plot([P7(1) P11(1) P10(1)],[P7(2) P11(2) P10(2)],'color',[244/255, 142/255 164/255],'LineWidth',2)
            plot([P_genou(1) P11(1)],[P_genou(2) P11(2)],'k:','LineWidth',4)
            hold off
            axis('equal')


            function plotcircle(C,r)
                % draw a circle of radius R around C
                a=0:0.01:2*pi;
                x=C(1)+ r*cos(a);
                y=C(2)+ r*sin(a);
                plot(x,y,'k','LineWidth',2);
            end

            function plotcirclepoint(C,r)
                % draw a circle of radius R around C
                a=0:0.01:2*pi;
                x=C(1)+ r*cos(a);
                y=C(2)+ r*sin(a);
                plot(x,y,'--','color','k');
            end
        end

        function plotcorps(obj,P_hanche,q)
            hauteur_corps = 4*obj.corps_joint_h;
            longueur_corps = 0.7*hauteur_corps;
            P6a=P_hanche-0.7*(hauteur_corps/2)*[-sin(q(1)+q(2)+q(3)+q(4));cos(q(1)+q(2)+q(3)+q(4))];
            P6c=P_hanche+0.1*(hauteur_corps/2)*[-sin(q(1)+q(2)+q(3)+q(4));cos(q(1)+q(2)+q(3)+q(4))];
            plot([P6a(1) P6c(1)],[P6a(2) P6c(2)],'color',[112/255 113/255 115/255] ,'LineWidth',10);

            %hip position
            a=0:0.01:2*pi;
            r=0.3;
            plot(P_hanche(1)+ r*cos(a),P_hanche(2)+ r*sin(a),'k','LineWidth',3); % Dessin Point de la hanche

            %body position
            Phi=-32.1*pi/180;
            V1=[-sin(q(1)+q(2)+q(3)+q(4)+Phi);cos(q(1)+q(2)+q(3)+q(4)+Phi)];
            V2=[cos(q(1)+q(2)+q(3)+q(4)+Phi);sin(q(1)+q(2)+q(3)+q(4)+Phi)];
            P6m=P_hanche+0.3*hauteur_corps*V1-0.2*longueur_corps*V2;
            P6b=P6m-0.5*(hauteur_corps)*V1;
            P8=P6m+0.5*(hauteur_corps)*V1;
            P9=P6m+0.8*(hauteur_corps)*V1+(longueur_corps)*V2;
            P10=P6m-0.8*(hauteur_corps)*V1+(longueur_corps)*V2;
            plot([P6b(1) P8(1) P9(1) P10(1) P6b(1)],[P6b(2) P8(2) P9(2) P10(2) P6b(2)],'color',[112/255 113/255 115/255] ,'LineWidth',2) %dessin du corps

        end

        function plotCoM(obj,q)
            P_CoM=obj.CoM(q);
            plot([P_CoM(1)-0.15 P_CoM(1)+0.15],[P_CoM(2), P_CoM(2)],'color',[255/255 164/255 36/255],'LineWidth',1), %jaune
            plot([P_CoM(1) P_CoM(1)],[P_CoM(2)-0.15, P_CoM(2)+0.15],'color',[255/255 164/255 36/255],'LineWidth',1), %jaune
            a=0:0.01:2*pi;
            r=0.10;
            plot(P_CoM(1)+ r*cos(a),P_CoM(2)+ r*sin(a),'color',[255/255 164/255 36/255],'LineWidth',1);
            r=0.2;
            plot(P_CoM(1)+ r*cos(a),P_CoM(2)+ r*sin(a),'color',[255/255 164/255 36/255],'LineWidth',3);
            plot([P_CoM(1) P_CoM(1)],[0 P_CoM(2)],'--','color',[255/255 164/255 36/255]), %jaune
        end


        function dessin_patte_point(obj,q)
            P3=[0;abs(obj.rP)];
            plot([obj.doigt_lx1 obj.doigt_lx2 P3(1) obj.doigt_lx1], [0 0 P3(2) 0],'k','LineWidth',2)
            hold on

            P_cheville= P3 + obj.L_tmt*[-sin(q(1)); cos(q(1))];
            P_genou= P_cheville + obj.L_tbt*[-sin(q(1)+q(2)); cos(q(1)+q(2))];
            P_hanche= P_genou + obj.L_f*[-sin(q(1)+q(2)+q(3)); cos(q(1)+q(2)+q(3))];
             plot([P3(1) P_cheville(1) P_genou(1) P_hanche(1)],[P3(2) P_cheville(2) P_genou(2) P_hanche(2)],':','Color','k','LineWidth',2)
            plotcircle(P_cheville,obj.rA)
            plotcircle(P_genou,obj.rK1)
            plotcircle(P3,obj.rP)
            plotcirclepoint(P_cheville,obj.rA/2)
            grid

            obj.plotcorpspoint(P_hanche,q)
            obj.plotCoM(q)

            %cable yellow
            if obj.rP>0
                plot([obj.rP obj.rP+obj.doigt_lx2/3],[obj.rP 2*obj.rP/3],'--','color',[255/255 164/255 36/255],'LineWidth',2)
                beta = 0;
            else
                plot([0 abs(obj.rP)+obj.doigt_lx2/3],[0 0],'--','color',[255/255 164/255 36/255],'LineWidth',2)
                beta=pi/2;
            end
            lo=sqrt((obj.L_tmt)^2-(obj.rP+obj.rA)^2);
            beta_2=-asin((obj.rA+obj.rP)/obj.L_tmt);
            R1=(obj.rA*P3+obj.rP*P_cheville)/(obj.rP+obj.rA)-obj.rP/(obj.rP+obj.rA)*[lo*sin(-q(1)+beta_2);lo*cos(-q(1)+beta_2)];
            R2=(obj.rA*P3+obj.rP*P_cheville)/(obj.rP+obj.rA)+obj.rA/(obj.rP+obj.rA)*[lo*sin(-q(1)+beta_2);lo*cos(-q(1)+beta_2)];
            plot([R1(1) R2(1)],[R1(2) R2(2)],'--','color',[255/255 164/255 36/255],'LineWidth',2)
            a=beta:(q(1)-beta_2-beta)/100:q(1)-beta_2;
            plot(P3(1)+obj.rP*cos(a),P3(2)+obj.rP*sin(a),'--','color',[255/255 164/255 36/255],'LineWidth',2); %jonction segment au niveau du pied
            beta=asin((obj.rK1)/(3/4*obj.L_tbt));
            R3=P_cheville-obj.rA*[cos(q(1)+q(2)-beta);sin(q(1)+q(2)-beta)];
            R4y=(P_cheville+3*P_genou)/4;
            plot([R3(1) R4y(1)],[R3(2) R4y(2)],'--','color',[255/255 164/255 36/255],'LineWidth',2)
            a=q(1)-beta_2:(q(2)-beta+beta_2)/100:q(1)+q(2)-beta;
            plot(P_cheville(1)- obj.rA*cos(a),P_cheville(2)- obj.rA*sin(a),'--','color',[255/255 164/255 36/255],'LineWidth',1);

            %green cable
            hv=0.9;
            beta=asin((obj.rA)/(hv*obj.L_tbt));
            R3=P_cheville-obj.rA*[cos(q(1)+q(2)-beta);sin(q(1)+q(2)-beta)];
            R4=((1-hv)*P_cheville+hv*P_genou);
            plot([R3(1) R4(1)],[R3(2) R4(2)],'--','color',[0 145/255 69/255],'LineWidth',2)
            hp=0.6;
            beta2=-asin((obj.rA)/(hp*obj.L_tmt));
            R3=P_cheville-obj.rA*[cos(q(1)-beta2);sin(q(1)-beta2)];
            R4=((1-hp)*P3+hp*P_cheville);
            plot([R3(1) R4(1)],[R3(2) R4(2)],'--','color',[0 145/255 69/255],'LineWidth',2)
            a=q(1)-beta2:(q(2)-beta+beta2)/100:q(1)+q(2)-beta;
            plot(P_cheville(1)- obj.rA*cos(a),P_cheville(2)- obj.rA*sin(a),'--','color',[0 145/255 69/255],'LineWidth',2);

            %blue cable
            hv=0.9;
            beta=asin((-obj.rA/2)/(hv*obj.L_tbt));
            R3=P_cheville+obj.rA/2*[cos(q(1)+q(2)-beta);sin(q(1)+q(2)-beta)];
            R4=((1-hv)*P_cheville+hv*P_genou);
            plot([R3(1) R4(1)],[R3(2) R4(2)],'--','color',[65/255, 105/255, 225/255],'LineWidth',2)
            hp=0.6;
            beta2=-asin((-obj.rA/2)/(hp*obj.L_tmt));
            R3=P_cheville+obj.rA/2*[cos(q(1)-beta2);sin(q(1)-beta2)];
            R4=((1-hp)*P3+hp*P_cheville);
            a=-pi()/2:pi()/180:q(1)+q(2)-beta;
            plot(P_cheville(1)+ obj.rA/2*cos(a),P_cheville(2)+ obj.rA/2*sin(a),'--','color',[65/255, 105/255, 225/255],'LineWidth',2);


            %red cable
            Dist=sqrt(obj.L_f^2+obj.corps_joint_h^2-2*obj.L_f*obj.corps_joint_h*cos(q(4)));
            alpha=q(1)+q(2)+q(3)+atan2(-obj.corps_joint_h*sin(q(4)),obj.L_f-obj.corps_joint_h*cos(q(4)));
            P7=P_genou+Dist*[-sin(alpha);cos(alpha)];
            a4=+atan2(-obj.corps_joint_h*sin(q(4)),obj.L_f-obj.corps_joint_h*cos(q(4)))-acos(obj.rK1/Dist);
            gamma2=pi/2+ q(1)+q(2)+q(3)+a4;
            a=q(1)+q(2):abs((gamma2-q(1)+q(2)))/100:gamma2;
            plot(P_genou(1)+ obj.rK1*cos(a),P_genou(2)+ obj.rK1*sin(a),'--','color',[1 0 46/255],'LineWidth',2);
            R8=P_genou+obj.rK1*[cos(gamma2);sin(gamma2)];
            plot([P7(1) R8(1)],[P7(2) R8(2)],'--','color',[1 0 46/255],'LineWidth',2)

            %pink cable
            P10=P_genou-obj.L_fibula*(P_cheville-P_genou)/obj.L_tmt;
            theta=fminbnd(@(theta) obj.lenght(q,theta),q(3),pi);
            P11=P_genou+obj.rK2*[-sin((q(1)+q(2)+theta));cos((q(1)+q(2)+theta))];
            plot([P7(1) P11(1) P10(1)],[P7(2) P11(2) P10(2)],'--','color',[244/255, 142/255 164/255],'LineWidth',2)
            plot([P_genou(1) P11(1)],[P_genou(2) P11(2)],'k:','LineWidth',4)
            hold off
            axis('equal')

            function plotcircle(C,r)
                % draw a circle of radius R around C
                a=0:0.01:2*pi;
                x=C(1)+ r*cos(a);
                y=C(2)+ r*sin(a);
                plot(x,y,'k','LineWidth',2);
            end

            function plotcirclepoint(C,r)
                % draw a circle of radius R around C
                a=0:0.01:2*pi;
                x=C(1)+ r*cos(a);
                y=C(2)+ r*sin(a);
                plot(x,y,'--','color','k');
            end
        end

        function plotcorpspoint(obj,P_hanche,q)
            hauteur_corps = 4*obj.corps_joint_h;
            longueur_corps = 0.7*hauteur_corps;
            P6a=P_hanche-0.7*(hauteur_corps/2)*[-sin(q(1)+q(2)+q(3)+q(4));cos(q(1)+q(2)+q(3)+q(4))];
            P6c=P_hanche+0.1*(hauteur_corps/2)*[-sin(q(1)+q(2)+q(3)+q(4));cos(q(1)+q(2)+q(3)+q(4))];
            plot([P6a(1) P6c(1)],[P6a(2) P6c(2)],':','color',[112/255 113/255 115/255] ,'LineWidth',10);

            %hip point
            a=0:0.01:2*pi;
            r=0.3;
            plot(P_hanche(1)+ r*cos(a),P_hanche(2)+ r*sin(a),'k','LineWidth',3); % Dessin Point de la hanche

            %body
            Phi=-32.1*pi/180;
            V1=[-sin(q(1)+q(2)+q(3)+q(4)+Phi);cos(q(1)+q(2)+q(3)+q(4)+Phi)];
            V2=[cos(q(1)+q(2)+q(3)+q(4)+Phi);sin(q(1)+q(2)+q(3)+q(4)+Phi)];
            P6m=P_hanche+0.3*hauteur_corps*V1-0.2*longueur_corps*V2;
            P6b=P6m-0.5*(hauteur_corps)*V1;
            P8=P6m+0.5*(hauteur_corps)*V1;
            P9=P6m+0.8*(hauteur_corps)*V1+(longueur_corps)*V2;
            P10=P6m-0.8*(hauteur_corps)*V1+(longueur_corps)*V2;
            plot([P6b(1) P8(1) P9(1) P10(1) P6b(1)],[P6b(2) P8(2) P9(2) P10(2) P6b(2)],'--','color',[112/255 113/255 115/255] ,'LineWidth',2) %dessin du corps

        end

    end
end