clear variables
numberOfVariables=3;
[R1,fval] = gamultiobj(@simple_multiobjective,numberOfVariables,[],[],[],[],[],[],@mycon);
disp(R1);
[m,n]=size(R1);
for i=1:m
    D=1; A=D./R1(i,2); C= D./R1(i,3); B= sqrt(A.^2+C.^2+D.^2-2*R1(i,1).*A.*C);
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
end
function tr=simple_multiobjective(R)
    tr(1)=abs((2/R(3)^2)+(1/R(2)^2)-(2*R(1)/(R(2)*R(3))));
    tr(2)=abs(1/R(2));
end
function [j,jeq]=mycon(x)
r=0.2;
    j(1)=abs(((1/x(3))-(x(1)/x(2))+(x(3)/x(2)))/sqrt(1+(1/x(2)^2)+(1/x(3)^2)-(2*x(1)/(x(2)*x(3)))))-sin(r);
    j(1)=abs(((1/x(3))-(x(1)/x(2))-(x(3)/x(2)))/sqrt(1+(1/x(2)^2)+(1/x(3)^2)-(2*x(1)/(x(2)*x(3)))))-sin(r);
    jeq=[];
end
% function tr=simple_multiobjective(R)
%     tr(1)=abs((2/R(3)^2)+1-(2*R(1)/(R(2)*R(3))));
%     tr(2)=abs(R(2)-0.01);
% end
% function [j,jeq]=mycon(x)    
%     j(1)=abs(((1/x(3))-(x(1)/x(2))+(x(3)/x(2)))/sqrt(1+(1/x(2)^2)+(1/x(3)^2)-(2*x(1)/(x(2)*x(3)))))-1;
%     j(2)=abs(((1/x(3))-(x(1)/x(2))-(x(3)/x(2)))/sqrt(1+(1/x(2)^2)+(1/x(3)^2)-(2*x(1)/(x(2)*x(3)))))-1;
%     j(3)=-(1/x(2)^2+1/x(3)^2+1-2*x(1)/(x(2)*x(3)));
%     jeq=[];
% end
