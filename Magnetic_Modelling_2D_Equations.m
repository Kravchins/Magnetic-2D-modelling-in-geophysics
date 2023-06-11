function[T,I,D,J] = Magnetic_Modelling_2D_Equations(Polygon,x,mu_0,Jr,Ir,Dr,Ji,Ii,Di,C)

%% This program is used to create the polygon and output the total magnetic field
%% Any use of this code MUST refer to the original publication:
% Kravchinsky, V.A., Hnatyshin, D., Lysak, B., Alemie, W., 2019. 
% Computation of magnetic anomalies caused by two‐dimensional structures of arbitrary shape: Derivation and Matlab implementation. 
% Geophysical Research Letters, 46(13), 7345-7351, https://doi.org/10.1029/2019GL082767

    del=acosd(cosd(Ii)*cosd(Di)*cosd(Ir)*cosd(Dr)+cosd(Ii)*sind(Di)*cosd(Ir)*sind(Dr)+sind(Ii)*sind(Ir));
    del=180-del;
    cos_del=cosd(del);
    
    J=sqrt((Ji.^2)+(Jr.^2)-(2*Ji*Jr*cos_del)); % Total magnetization (induced + remnant) already converted from A/m to nT
    I=asind((Ji*sind(Ii)+Jr*sind(Ir))/J);      % Inclination for total magnetization

    JH=J*cosd(I);  % Horizontal component of total magnetization
 
    if Ji*cosd(Ii)*cosd(Di)+Jr*cosd(Ir)*cosd(Dr) == 0
      D = 0;
    else
      D = acosd((Ji*cosd(Ii)*cosd(Di)+Jr*cosd(Ir)*cosd(Dr))/JH); % Declination for total magnetization
                 
    end
    
    I=real(I);
    C=real(C);
    D=real(D); 

polylength = length(Polygon(:,1)); 

Jx = J*cosd(I)*cosd(C-D); % X component it is also equal to J*cos(I)*sin(D) since C=90
Jz = J*sind(I);           % Z component
T = [];

for n = 1:length(x)
    for i=1:polylength
        if i == polylength
            x1 = Polygon(i,1)-x(n);
            x2 = Polygon(1,1)-x(n);
            z1 = Polygon(i,2);
            z2 = Polygon(1,2);   
        else
            x1 = Polygon(i,1)-x(n);
            x2 = Polygon(i+1,1)-x(n);
            z1 = Polygon(i,2);
            z2 = Polygon(i+1,2);
        end
       
        if z2-z1 == 0
            V(n,i) = 0;
            H(n,i) = 0;
        else 
            g(n,i) = (x2-x1)/(z2-z1);
            if (x1-z1*(g(n,i))) >= 0
              d(n,i) = 1;
            else
              d(n,i) = -1;
            end
            
          alfa1(n,i) = atan2((z1+x1*g(n,i)),abs(x1-z1*g(n,i)));
          alfa2(n,i) = atan2((z2+x2*g(n,i)),abs(x2-z2*g(n,i)));
        
        r1(n,i) = sqrt((x1)^2+z1^2);
        r2(n,i) = sqrt((x2)^2+z2^2);
        x21=x2-x1;
        z21=z2-z1;
        gamz(n,i)=z21/sqrt(x21.^2+z21.^2);
        gamx(n,i)=x21/sqrt(x21.^2+z21.^2);
         
        s=-1;   %Introducing a negative sign to use the clockwise order of points for the polygon
        Q(n,i) = s*((gamz(n,i)^2*log(r2(n,i)/r1(n,i))-(d(n,i)*gamz(n,i)*gamx(n,i)*(alfa2(n,i)-alfa1(n,i)))));    
        P(n,i) = s*((gamz(n,i)*gamx(n,i)*log(r2(n,i)/r1(n,i)) + d(n,i)*gamz(n,i)^2*(alfa2(n,i)-alfa1(n,i))));
%       It can be also written as
%       P(n,i) = s*((z2-z1)^2*d(n,i)*(alfa2(n,i)-alfa1(n,i))/((z2-z1)^2+(x2-x1)^2) + (x2-x1)*(z2-z1)*log(r2(n,i)/r1(n,i))/((z2-z1)^2+(x2-x1)^2));
%       Q(n,i) = s*((z2-z1)^2*log(r2(n,i)/r1(n,i))/((z2-z1)^2+(x2-x1)^2)-(z2-z1)*(x2-x1)*d(n,i)*(alfa2(n,i)-alfa1(n,i))/((z2-z1)^2+(x2-x1)^2));    
            
        V(n,i) = 2*(Jx*Q(n,i) - Jz*(P(n,i)))/(4*pi);
        H(n,i) = 2*(Jx*P(n,i) + Jz*(Q(n,i)))/(4*pi);
            
        end    
end
    T(n) = (sum(V(n,:))*sind(I) + sum(H(n,:))*cosd(I)*cosd(C-D));
    
end
display 'Any use of this code MUST refer to the original publication:'
display 'Kravchinsky, V.A., Hnatyshin, D., Lysak, B., Alemie, W., 2019.'
display 'Computation of magnetic anomalies caused by two‐dimensional structures of arbitrary shape: Derivation and Matlab implementation.'
display 'Geophysical Research Letters, 46(13), 7345-7351, https://doi.org/10.1029/2019GL082767'
end
