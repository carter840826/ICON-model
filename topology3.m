function [t1,t2,t3] = topology3( M,r,p,q,x,x1,x2,x3,z1,z2,z3 )
%UNTITLED Summary of this function goes here
%t = (0:0.008:12);
t = 1:M;

for i = 1:M-1
%y1(i) = (x1(i+1)-x1(i))/0.008;
%y2(i) = (x2(i+1)-x2(i))/0.008;
%y3(i) = (x3(i+1)-x3(i))/0.008;
y(i) = (x(i+1)-x(i));
y1(i) = (x1(i+1)-x1(i));
y2(i) = (x2(i+1)-x2(i));
y3(i) = (x3(i+1)-x3(i));

end
t1(1)=x1(1);t2(1)=x2(1);t3(1)=x3(1);
for m = 1:M-1
    K1=[];
        %for K1
    for k = 1:r
            A1 = [sin(k*(t1(m))),cos(k*(t1(m)))];
            K1 = [K1 A1];
    end
    for k = 1:r
            A1 = [sin(k*x(m))*sin(k*t1(m)),sin(k*t2(m))*sin(k*t1(m)),sin(k*x3(m))*sin(k*t1(m))];
            %A1 = [sin(k*(t2(m)-t1(m))),sin(k*(t3(m)-t1(m)))];
            K1 = [K1 A1];
    end
    for k = 1:r
            A1 = [cos(k*x(m))*cos(k*t1(m)),cos(k*t2(m))*cos(k*t1(m)),cos(k*x3(m))*cos(k*t1(m))];
            %A1 = [cos(k*(t2(m)-t1(m))),cos(k*(t3(m)-t1(m)))];
            K1 = [K1 A1];
    end
    %for K2
    K2=[];
    for k = 1:p
            A2 = [sin(k*(t2(m))),cos(k*(t2(m)))];
            K2 = [K2 A2];
    end
    for k = 1:p
            A2 = [sin(k*x(m))*sin(k*t2(m)),sin(k*t1(m))*sin(k*t2(m)),sin(k*t3(m))*sin(k*t2(m))];
            %A2 = [sin(k*(t1(m)-t2(m))),sin(k*(t3(m)-t2(m)))];
            K2 = [K2 A2];
    end
    for k = 1:p
            A2 = [cos(k*x(m))*cos(k*t2(m)),cos(k*t1(m))*cos(k*t2(m)),cos(k*t3(m))*cos(k*t2(m))];
            %A2 = [cos(k*(t1(m)-t2(m))),cos(k*(t3(m)-t2(m)))];
            K2 = [K2 A2];
    end
        %for K3
    K3=[1];
    %for k = 1:q
            %213
            %A3 = [sin(k*(x3(m))),cos(k*(x3(m)))];
            %A1 = [sin(k*(x1(m)-x1(m))),sin(k*(x2(m)-x1(m))),sin(k*(x3(m)-x1(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            %K3 = [K3 A3];
   % end
    for k = 1:q
        A3 = [sin(k*(t1(m)-t3(m))),sin(k*(t2(m)-t3(m)))];    
        %A3 = [sin(k*x1(m))*sin(k*x3(m)),sin(k*x2(m))*sin(k*x3(m))];
            %A3 = [sin(k*(x1(m)-x3(m))),sin(k*(x2(m)-x3(m))),sin(k*(x3(m)-x3(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
        K3 = [K3 A3];
    end
    for k = 1:q
        A3 = [cos(k*(t1(m)-t3(m))),cos(k*(t2(m)-t3(m)))];    
        %A3 = [cos(k*x1(m))*cos(k*x3(m)),cos(k*x2(m))*cos(k*x3(m))];
            %A3 = [cos(k*(x1(m)-x3(m))),cos(k*(x2(m)-x3(m))),cos(k*(x3(m)-x3(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K3 = [K3 A3];
    end
        
        d1(m) = K1*z1;
        d2(m) = K2*z2;
        d3(m) = K3*z3;
        t1(m+1) = t1(m)+d1(m);
        t2(m+1) = t2(m)+d2(m);
        t3(m+1) = t3(m)+d3(m);
end


subplot(4,1,1), plot(t,x(1:M));ylim([10 35]);title('天氣預報');
subplot(4,1,2), ax=plotyy(t,x1(1:M),t,t1, 'plot', 'plot');ylim(ax(1),[10 20]);ylim(ax(2),[10 20]);title('土壤溫度');
subplot(4,1,3), bx=plotyy(t,x2(1:M),t,t2, 'plot', 'plot');ylim(bx(1),[10 20]);ylim(bx(2),[10 20]);title('室內溫度');
subplot(4,1,4), cx=plotyy(t,x3(1:M),t,t3, 'plot', 'plot');ylim(cx(1),[0.2 0.35]);ylim(cx(2),[0.2 0.35]);title('土壤含水量');
end