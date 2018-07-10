clc;
clear variables
r1=input('enter range: 0 to ');
t=ceil(max(3*r1/(2*pi),3));
format='Please select number of presicion points more than %d \n';
fprintf(format,t);
prompt='enter number of presicion points ';
n=input(prompt);
n=ceil(n);
if(n>t)
f=r1;
theta_s= 0;  
theta_f= f;
phi_s= 0;
phi_f= f;
i=1:n;
theta = theta_s + (theta_f - theta_s)*(i-1)/(n-1); 
x_s = 0;
x_f = f;
r_x=(theta_f-theta_s)/(x_f-x_s);
x = x_s + (theta-theta_s)/r_x;
y= cos(x);
t0=0;
s=zeros(100,1);
for r_y=100:-1:1
phi=phi_s+ r_y*(y- y(1));
    phi=deg2rad(phi);
    sum1=sum(cos(theta));
    sum2=sum(cos(phi));
    sum3=sum(cos(phi).^2);
    sum4=sum(cos(theta).^2);
    sum5=sum(cos(theta-phi));
    sum6=sum(cos(theta).*cos(phi));
    sum7=sum(cos(theta-phi).*cos(phi));
    sum8=sum(cos(theta-phi).*cos(theta));
    A1=[1,sum2,-sum1;sum2,sum3,-sum6;sum1,sum6,-sum4];
    B1=[sum5,sum7,sum8];
    R= A1\B1';
    
    D=5; A=D/R(2); C= D/R(3); B= sqrt(A^2+C^2+D^2-2*R(1)*A*C);
    t1=B+D-A-C;
    t2=D+C-A-B;
    t3=B+C-A-D;
    A=A/4;B=B/4;C=C/4;D=D/4;
    E = sqrt(A^2 + D^2 - 2*A*D*cos(theta));
    j=(E.^2 + C^2 - B^2)./(2*E*C);
    if(abs(max(j))<1)
        %disp('hi');
        if(abs((C^2 + B^2 - E.^2)./(2*B*C))<1)
        %disp('not complex');
    alfa = asin(A*sin(theta)./E);
    beta = acos((E.^2 + C^2 - B^2)./(2*E*C));
    tran=acos((C^2 + B^2 - E.^2)./(2*B*C)); 
    t=[-abs(t1) -abs(t2) abs(t3);abs(t1) abs(t2) abs(t3);abs(t1) -abs(t2) -abs(t3);-abs(t1) abs(t2) -abs(t3)]';
    if(t1==t(1,1) && t2== t(2,1) && t3== t(3,1))
        %disp('crank-crank');
    elseif(t1==t(1,2) && t2== t(2,2) && t3== t(3,2))
         %disp('crank-rocker');
    elseif(t1==t(1,3) && t2== t(2,3) && t3== t(3,3))
         %disp('rocker-crank');
    elseif(t1==t(1,4) && t2== t(2,4) && t3== t(3,4))
         %disp('rocker-rocker');
    else
        %disp('non-grashof');
        s(r_y)=n;
%         disp(r_y);
    end
    else
        %disp('complex');
        t0(r_y)=1;
        %disp(j);
        continue;
        end
    else
        %disp('complex');
        t0(r_y)=1;
        %disp(j);
        continue;
    end
end
[m,n]=size(s);
for i=1:n
    if(s(i)==0)
            r_y=i;
        
    end
end
phi=phi_s+ r_y*(y- y(1));
    phi=deg2rad(phi);
    sum1=sum(cos(theta));
    sum2=sum(cos(phi));
    sum3=sum(cos(phi).^2);
    sum4=sum(cos(theta).^2);
    sum5=sum(cos(theta-phi));
    sum6=sum(cos(theta).*cos(phi));
    sum7=sum(cos(theta-phi).*cos(phi));
    sum8=sum(cos(theta-phi).*cos(theta));
    A1=[1,sum2,-sum1;sum2,sum3,-sum6;sum1,sum6,-sum4];
    B1=[sum5,sum7,sum8];
    R= A1\B1';
 R=[   -0.5708   -0.0828    0.1477];

% R=[  -0.6394   -0.2253    0.3776];
    R=R';

    D=1; A=D/R(2); C= D/R(3); B= sqrt(A^2+C^2+D^2-2*R(1)*A*C);
%     disp(B^2+C^2-D^2);
[m,n]=size(x);
xo=linspace(0,f,n);
yo=cos(xo);

disp([A B C D]);
    t1=B+D-A-C;
    t2=D+C-A-B;
    t3=B+C-A-D;
    t=[-abs(t1) -abs(t2) abs(t3);abs(t1) abs(t2) abs(t3);abs(t1) -abs(t2) -abs(t3);-abs(t1) abs(t2) -abs(t3)]';
    if(t1==t(1,1) && t2== t(2,1) && t3== t(3,1))
        disp('crank-crank');
    elseif(t1==t(1,2) && t2== t(2,2) && t3== t(3,2))
         disp('crank-rocker');
    elseif(t1==t(1,3) && t2== t(2,3) && t3== t(3,3))
         disp('rocker-crank');
    elseif(t1==t(1,4) && t2== t(2,4) && t3== t(3,4))
         disp('rocker-rocker');
    else
        disp('Non-grashof');
    end
t=theta;
xx = 0:.1:f;
P1 = [0;0];
P4 = D*[1;0];
%disp((C.^2 + B.^2 - E.^2)./(2*B*C));
P2 = A*[cos(theta); sin(theta)]; 
E = sqrt(A^2 + D^2 - 2*A*D*cos(theta));
alfa = asin(A*sin(theta)./E);
beta = acos((E.^2 + C^2 - B^2)./(2*E*C));
tran=acos((C.^2 + B.^2 - E.^2)./(2*B*C));
P3 = [D - C*cos(alfa+beta); C*sin(alfa+beta)];
disp(tran);

P3_x = P3(1,:);
P3_y = P3(2,:);
P3_vx = diff(P3_x)./diff(t);
P3_vy = diff(P3_y)./diff(t);
t1=t(1:n-1);
P3_ax = diff(P3_vx)./diff(t1);
P3_ay = diff(P3_vy)./diff(t1);

P3_v = sqrt(P3_vx.^2 + P3_vy.^2);
P3_a = sqrt(P3_ax.^2 + P3_ay.^2);

yy = spline(x,y,xx);
y1=y-spline(x,y,x);

for i=1:length(t)-1
   subplot(2,2,3);
   set(gca,'XLim',[0 f],'YLim',[-1 1]);
   plot(xo,yo,x(1:i),y(1:i),'o',xx,yy,x(1:i),y1(1:i),'o');
   grid on
   grid minor
   hold on
   set(gca,'XLim',[0 f],'YLim',[-inf inf]);
   plot(theta(1:i),tran(1:i),'r');
   grid on
   grid minor
   hold on
   ani = subplot(2,2,1);
   P1_circle = viscircles(P1',0.05);
   P2_circle = viscircles(P2(:,i)',0.05);
   P3_circle = viscircles(P3(:,i)',0.05);
   P4_circle = viscircles(P4',0.05); 
   
   A_bar = line([P1(1) P2(1,i)],[P1(2) P2(2,i)]);
   B_bar = line([P2(1,i) P3(1,i)],[P2(2,i) P3(2,i)]);
   C_bar = line([P3(1,i) P4(1)],[P3(2,i) P4(2)]);
   
   axis(ani,'equal');
   set(gca,'XLim',[-5 20],'YLim',[-5 20]);
   
   str1 = 'P3';
   str3='o';
   str4='.';
   str2 = ['Time elapsed: '  num2str(t(i)) ' s'];
   P3_text = text(P3(1,i),P3(2,i)+0.6,str1);
   P3_text1 = text(P3(1,i),P3(2,i),str3,'color','red');
   P1_text1 = text(P2(1,i),P2(2,i),str4,'color','red');
   Time = text(-2,6,str2);
   pause(0.005);
   if i<length(t1)
    delete(P1_circle);
    delete(P2_circle); 
    delete(P3_circle);
    delete(P4_circle);
    delete(A_bar);
    delete(B_bar);
    delete(C_bar);
    delete(P3_text);
    delete(Time);
    vel = subplot(2,2,2);
    plot(vel,t1(1:i),P3_v(1:i));
    set(vel,'XLim',[0 t(n-1)],'YLim',[0 inf]);
    xlabel(vel, 'Time (s)');
    ylabel(vel, 'Amplitude (m/s)');
    title(vel,'Speed of P3');
    grid on;
    acc = subplot(2,2,4);
    plot(acc,t1(1:i),P3_a(1:i));
    set(acc,'XLim',[0 t(n-1)],'YLim',[0 inf]);
    xlabel(acc, 'Time (s)');
    ylabel(acc, 'Amplitude (m/s^2)');
    title(acc,'Acceleration of P3');
    grid on;
   end
end
else
    disp('too less presecion points for this range');
end
