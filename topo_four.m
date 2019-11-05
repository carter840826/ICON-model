function [ t1,t2,t3,z1,z2,z3 ] = topo_four( M,r,p,q,x,x1,x2,x3 )
%UNTITLED2 Summary of this function goes here
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

%y1场だ
G1 = [];
for m = 1:M-1
        K=[];
        for k = 1:r
            %260
            A1 = [sin(k*(x1(m))),cos(k*(x1(m)))];
            %A1 = [sin(k*(x1(m)-x1(m))),sin(k*(x2(m)-x1(m))),sin(k*(x3(m)-x1(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A1];
        end
        for k = 1:r
            A1 = [sin(k*x(m))*sin(k*x1(m)),sin(k*x2(m))*sin(k*x1(m)),sin(k*x3(m))*sin(k*x1(m))];
            %A1 = [sin(k*(x1(m)-x1(m))),sin(k*(x2(m)-x1(m))),sin(k*(x3(m)-x1(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A1];
        end
        for k = 1:r
            A1 = [cos(k*x(m))*cos(k*x1(m)),cos(k*x2(m))*cos(k*x1(m)),cos(k*x3(m))*cos(k*x1(m))];
            %A1 = [cos(k*(x1(m)-x1(m))),cos(k*(x2(m)-x1(m))),cos(k*(x3(m)-x1(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A1];
        end
        G1 = [G1;K];
end
%よ猭
%z = inv(G'*G)*G'*s1';
%よ猭
%[U, S, V] = svd(G);
%AH = V*S'*U';
%Z = AH*s1';
%よ猭
[U1, S1, V1] = svd(G1'*G1);
for i = 1:length(S1(1,:))
    if S1(i,i) < 10^-15
        S1(i,i) = 0;
    else
        S1(i,i) = S1(i,i);
    end
end
z1 = inv(U1*S1*V1')*G1'*y1';
d1 = G1 *z1;
t1(1) = x1(1);
for i = 2:M
    t1(i) = t1(i-1) + d1(i-1);
    %t1(i) = t1(i-1) + d1(i-1)*0.008;
end

%y2场だ
G2 = [];
for m = 1:M-1
        K=[];
        for k = 1:p
            %213
            A2 = [sin(k*(x2(m))),cos(k*(x2(m)))];
            %A1 = [sin(k*(x1(m)-x1(m))),sin(k*(x2(m)-x1(m))),sin(k*(x3(m)-x1(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A2];
        end
        for k = 1:p
            A2 = [sin(k*x(m))*sin(k*x2(m)),sin(k*x1(m))*sin(k*x2(m)),sin(k*x3(m))*sin(k*x2(m))];
            %A2 = [sin(k*(x1(m)-x2(m))),sin(k*(x2(m)-x2(m))),sin(k*(x3(m)-x2(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A2];
        end
        for k = 1:p
            A2 = [cos(k*x(m))*cos(k*x2(m)),cos(k*x1(m))*cos(k*x2(m)),cos(k*x3(m))*cos(k*x2(m))];
            %A2 = [cos(k*(x1(m)-x2(m))),cos(k*(x2(m)-x2(m))),cos(k*(x3(m)-x2(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A2];
        end
        G2 = [G2;K];
end
[U2, S2, V2] = svd(G2'*G2);
for i = 1:length(S2(1,:))
    if S2(i,i) < 10^-15
        S2(i,i) = 0;
    else
        S2(i,i) = S2(i,i);
    end
end
z2 = inv(U2*S2*V2')*G2'*y2';
d2 = G2 *z2;
t2(1) = x2(1);
for i = 2:M
    t2(i) = t2(i-1) + d2(i-1);
    %t2(i) = t2(i-1) + d2(i-1)*0.008;
end

%y3场だ
G3 = [];
for m = 1:M-1
        K=[1];
        %for k = 1:q
            %213
            %A3 = [sin(k*(x3(m))),cos(k*(x3(m)))];
            %A1 = [sin(k*(x1(m)-x1(m))),sin(k*(x2(m)-x1(m))),sin(k*(x3(m)-x1(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            %K = [K A3];
        %end
        for k = 1:q
            A3 = [sin(k*(x1(m)-x3(m))),sin(k*(x2(m)-x3(m)))];
            %A3 = [sin(k*x1(m))*sin(k*x3(m)),sin(k*x2(m))*sin(k*x3(m))];
            %A3 = [sin(k*(x1(m)-x3(m))),sin(k*(x2(m)-x3(m))),sin(k*(x3(m)-x3(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A3];
        end
        for k = 1:q
            A3 = [cos(k*(x1(m)-x3(m))),cos(k*(x2(m)-x3(m)))];
            %A3 = [cos(k*x1(m))*cos(k*x3(m)),cos(k*x2(m))*cos(k*x3(m))];
            %A3 = [cos(k*(x1(m)-x3(m))),cos(k*(x2(m)-x3(m))),cos(k*(x3(m)-x3(m)))];
            %P = [legendreP(p,avet3(m))*legendreP(r,avet3(m)),legendreP(p,avet3(m))*legendreP(r,avet5(m)),legendreP(p,avet3(m))*legendreP(r,avet6(m))];
            K = [K A3];
        end
        G3 = [G3;K];
end
[U3, S3, V3] = svd(G3'*G3);
for i = 1:length(S3(1,:))
    if S3(i,i) < 10^-15
        S3(i,i) = 0;
    else
        S3(i,i) = S3(i,i);
    end
end
z3 = inv(U3*S3*V3')*G3'*y3';
d3 = G3 *z3;
t3(1) = x3(1);
for i = 2:M
    t3(i) = t3(i-1) + d3(i-1);
    %t3(i) = t3(i-1) + d3(i-1)*0.008;
end


subplot(4,1,1), plot(t,x(1:M));ylim([10 35]);title('ぱ箇厨');
subplot(4,1,2), ax=plotyy(t,x1(1:M),t,t1, 'plot', 'plot');ylim(ax(1),[10 35]);ylim(ax(2),[10 35]);title('腫放');
subplot(4,1,3), bx=plotyy(t,x2(1:M),t,t2, 'plot', 'plot');ylim(bx(1),[10 35]);ylim(bx(2),[10 35]);title('ず放');
subplot(4,1,4), cx=plotyy(t,x3(1:M),t,t3, 'plot', 'plot');ylim(cx(1),[0 0.5]);ylim(cx(2),[0 0.5]);title('腫秖');
end



