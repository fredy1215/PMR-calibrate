clc,clear,close all
%%
t = 0:15000;
M0 = 1;
TE = [60,30,50,25,40,0];
T2 = 195; % [ms]
n = 208;

TR = 11;
T1 = 1895; %1895
M = [];
for q = 1:numel(TE)
M1{q} = M0*exp(-t(1:TE(q)+1)/T2); % t=1:61
M_1{q} = -M1{q}(TE(q)+1); % t = 62

M(q,1,1,1:11) = M_1{q}*exp(-(1:11)/T1)+M0*(1-exp(-(1:11)/T1)); %t =63:72
M(q,1,2,1) = M(q,1,1,11)*cosd(8);
for c = 1:n-1
    for cc = 1:10
    M(q,1,c+1,cc+1)=M(q,1,c+1,1)*exp(-(cc+1)/T1)+M0*(1-exp(-(cc+1)/T1));
    end
    if c<n-1
    M(q,1,c+2,1) = M(q,1,c+1,11)*cosd(8);
    end
end
for k = 1:4
    M(q,k+1,1,1) = -M(q,k,208,11)*cosd(8);
    for c = 1:n
        for cc = 1:10
            M(q,k+1,c,cc+1) = M(q,k+1,c,1)*exp(-(cc+1)/T1)+M0*(1-exp(-(cc+1)/T1));
        end
        if c <n
            M(q,k+1,c+1,1) = M(q,k+1,c,11)*cosd(8);
        end
    end
end
end
 u = 1;
for q = 1:numel(TE)
    for k = 1:5
        for c = 1:n
    Signal(q,u) = M(q,k,c,11)*sind(8);
    u = u+1;
        end
    end
    u = 1;
end
            
%% continues signal

u = 1;
Mtot = [];
for q = 1:numel(TE)
for k = 1:5
for c = 1:n
    for cc = 1:11
Mtot(u) = M(q,k,c,cc);
u = u+1;
    end
end
end
figure('name',sprintf('TE = %i',TE(q))) 
plot(t(1:TE(q)+1),M1{q})
hold on
plot(t(TE(q)+2),M_1{q})
plot(t(TE(q)+3:numel(Mtot)+TE(q)+2), Mtot)


hold off
u = 1;
Mtot = [];
end
%% messured signal
a = 0:11:11*208*5-1;
for q = 1:numel(TE)
figure('name',sprintf('signal TE = %i',TE(q)))
plot(a,Signal(q,:))
end





