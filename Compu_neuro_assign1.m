%% Question 1: 

clearvars;
close all;
r = 0.6;
T = 0.1;
K = 0.02;
I = 0.03;
logthr=@(x) -r*x.*(1-x/T).*(1-x/K);
Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;
x = zeros(1,length(t));
x(1) = 0.2;
for j=1:length(t)-1
k1x = logthr(x(j))+I;
ax = x(j)+k1x*dt;
k2x = logthr(ax)+I;
x(j+1)=x(j)+(k1x+k2x)*dt/2;
end
xx = -1:0.01:2;
figure(1)
hFig = figure(1);
set(hFig, 'Position', [40 400 1000 500]); subplot(1,2,1)
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
xlabel('t');
ylabel('V');
subplot(1,2,2)
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]); plot(xx,logthr(xx)+I,'linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',24);
xlabel('V');
ylabel('F');
MVN = 0;
if MVN == 1
xx = -1:0.01:2;
figure(101)
k = 1;
hFig = figure(101);
set(hFig, 'Position', [20 100 700 500])
for j=1:floor(length(t)/k)-1
subplot(1,2,1)
plot(t,x,'-b','linewidth',2);
hold on
plot(t(k*j),x(k*j),'or','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
xlabel('t');
ylabel('y');
pause(0.0001);
hold off
subplot(1,2,2)
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]); hold on
plot(xx,logthr(xx)+I,'linewidth',2); plot(x(k*j),0,'or','linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('X');
pause(0.0001);
hold off
end
end

%% Question 2 
%Iapp = 0.05
clearvars;
close all;
r = 1;
T = 0.25;
K = 1;
I = 0.05;
logthr=@(x) -r*x.*(1-x/T).*(1-x/K);
Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;
x = zeros(1,length(t));
x(1) = 0.01;
for j=1:length(t)-1
k1x = logthr(x(j))+I;
ax = x(j)+k1x*dt;
k2x = logthr(ax)+I;
x(j+1)=x(j)+(k1x+k2x)*dt/2;
end
xx = -1:0.01:2;
figure(1)
hFig = figure(1);
set(hFig, 'Position', [40 400 1000 500]); subplot(1,2,1)
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
xlabel('t');
ylabel('V');
subplot(1,2,2)
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]); plot(xx,logthr(xx)+I,'linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('V');
ylabel('F');
MVN = 0;
if MVN == 1
xx = -1:0.01:2;
figure(101)
k = 1;
hFig = figure(101);
set(hFig, 'Position', [20 100 700 500])
for j=1:floor(length(t)/k)-1 subplot(1,2,1);
plot(t,x,'-b','linewidth',2);
hold on plot(t(k*j),x(k*j),'or','linewidth',2); 
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',24);
xlabel('t');
ylabel('V');
pause(0.0001);
hold off
subplot(1,2,2)
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]); hold on
plot(xx,logthr(xx)+I,'linewidth',2);
plot(x(k*j),0,'or','linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',24);
xlabel('V');
pause(0.0001);
hold off
end
end

%% Question 2
% Iapp = 0

clearvars;
close all;
% ODE
% X' = -r (1-X/T) ( 1 - X/K) X + I;
r = 1;
T = 0.25;
K = 1;
I = 0;
logthr=@(x) -r*x.*(1-x/T).*(1-x/K);
Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;
x = zeros(1,length(t));
x(1) = 0.01;
for j=1:length(t)-1
k1x = logthr(x(j))+I;
ax = x(j)+k1x*dt;
k2x = logthr(ax)+I;
x(j+1)=x(j)+(k1x+k2x)*dt/2;
end
xx = -1:0.01:2;
figure(1)
hFig = figure(1);
set(hFig, 'Position', [40 400 1000 500]); subplot(1,2,1)
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
xlabel('t');
ylabel('V');
subplot(1,2,2)
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]); plot(xx,logthr(xx)+I,'linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('V');
ylabel('F');
MVN = 0;
if MVN == 1
xx = -1:0.01:2;
figure(101)
k = 1;
hFig = figure(101);
set(hFig, 'Position', [20 100 700 500])
for j=1:floor(length(t)/k)-1
subplot(1,2,1)
plot(t,x,'-b','linewidth',2);
hold on
plot(t(k*j),x(k*j),'or','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',24);
xlabel('t');
ylabel('V');
pause(0.0001);
hold off
subplot(1,2,2)
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]); hold on
plot(xx,logthr(xx)+I,'linewidth',2); plot(x(k*j),0,'or','linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',24);
xlabel('V');
pause(0.0001);
hold off
end
end

%% Question 4: 
clearvars;
close all;
r=1;
T=0.4;
K=1;
I_app=0.0;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.02;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.04;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.06;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.08;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.10;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.12;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app);
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.14
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.16
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.18;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.20;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold off
figure(1);
plot (I,V,'r','LineWidth',1.5);
xlabel('I_{app}');
ylabel('V_{eq}');
axis([0.2 4.5 -0.5 2.0])
yline(0,'--')

hold on
I_app=0.18;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.16;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.14;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.12;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.10;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.08;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.06;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.04;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.02;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold on
I_app=0.0;
I_app=[0,0.2];
DI=0:0.02:0.2;
X_in=0.25;
sort(I_app)
V_in=0;
f = @(I,V) -r * V *(1-(V/T)) * (1-(V/K)) + I;
tspan = [0 20];
[I,V] = ode45(f , tspan , V_in);
hold off
figure(2)
plot (I,V,'b','LineWidth',1);
xlabel('I_{app}');
ylabel('V_{eq}');
axis([0.2 4.5 -0.5 2.0])
yline(0,'--')
%% Question 5:
clear
close all
r = 1;
T = 0.4;
K = 1;
I = 0 : 0.02 : 0.2;
I_descendend = 0.2 : -0.02 : 0;
Veq = zeros(1, length(I));
VeqDesc = zeros(1, length(I));
logthr=@(x) -r*x.*(1-x/T).*(1-x/K);
x = 0;
dt = 0.01;
for i = 1 : length(I)
    total = 1;
    while total > 1e-6
        k1x = logthr(x) + I(i);
        ax = x + k1x*dt;
        k2x = logthr(ax) + I(i);
        x_new = x + (k1x+k2x)*dt/2;
        total = abs(x_new - x);
x = x_new; end
    Veq(i) = x;
end
for i = 1 : length(I_descendend)
    total = 1;
    while total > 1e-6
        k1x = logthr(x) + I_descendend(i);
        ax = x + k1x*dt;
        k2x = logthr(ax) + I_descendend(i);
        x_new = x + (k1x+k2x)*dt/2;
        total = abs(x_new-x);
        x = x_new;
end
    VeqDesc(i) = x_new;
end
figure(1)
hold on
plot(I, Veq,'G','linewidth',2);
plot(I_descendend, VeqDesc, 'R','linewidth',2);
hold off
xlabel('Iapp');
ylabel('Veq');
legend({'ascending', 'descending'})
