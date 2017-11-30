clear;
clc;
%%
Basics();
%% Part 1
tau_ir = 5000;
[Spikes,L] =  GenSpike();
[Spikes,V,g,xarr] = verNoP(Spikes,L,tau_ir);
plotspikes(Spikes)
%% Part 2
PSTH5 = PSTH_50itr(L,tau_ir);
suptitle('\tau_i_r = 5000ms')
%% Part 3
tau_ir = 1000;
PSTH1 = PSTH_50itr(L,tau_ir);
suptitle('\tau_i_r = 1000ms')
tau_ir = 3000;
PSTH3 = PSTH_50itr(L,tau_ir);
suptitle('\tau_i_r = 3000ms')
tau_ir = 10000;
PSTH10 = PSTH_50itr(L,tau_ir);
suptitle('\tau_i_r = 10000ms')
%% Functions
function Basics()
n = 700;
lamda  = 40+rand(1,1000)*(50-40);
p = lamda./1000;         %probabilty of spiking per 1 ms
Spikes = zeros(n,1000);
for i = 1:n
    for j = 1:1000
        if(rand() < p(j))
            Spikes(i,j) = 1;
        end
    end
end
SpikeTime = cell(320,1);
for i = 1:n
    SpikeTime{i,1} = find(Spikes(i,:));
end
%%
a = [10 20 40 80 160 320 700];
mse = zeros(1,length(a));
PSTH = zeros(length(a),1000);
for i = 1:length(a)
    PSTH(i,:) = 1000*mean(Spikes(1:a(i),:),1);
    mse(i) = sqrt(immse(lamda,PSTH(i,:)));
end
plot(a,mse,'x-');
ylabel('Mean Square Error');
xlabel('# of Repetitions');
%figure;
end
function PSTH1 = PSTH_50itr(L,tau_ir)
PSTH.sp = zeros(1,L);
PSTH.l4 = zeros(1,L);
for i = 1:50
    [Spikes,L] =  GenSpike();
    [Spikes,~,~,~] = verNoP(Spikes,L,tau_ir);
    PSTH.sp = PSTH.sp + Spikes.sp;
    PSTH.l4 = PSTH.l4 + Spikes.l4;
end
PSTH.sp = PSTH.sp * 1000/50;
PSTH.l4 = PSTH.l4 * 1000/50;
for i = 1:L/10
    PSTH1.sp(i) = mean(PSTH.sp(10*(i-1)+1:10*i));
    PSTH1.l4(i) = mean(PSTH.l4(10*(i-1)+1:10*i));
end
figure;
subplot(2,1,1)
plot(PSTH1.sp)
xlabel('t(x10ms)');
ylabel('PSTH_S_P');
subplot(2,1,2)
plot(PSTH1.l4)
xlabel('t(x10ms)');
ylabel('PSTH_L_4');
end
function [Spikes,V,g,xarr] = verNoP(Spikes,L,tau_ir)
Spikes.sp = zeros(1,L);
Spikes.l4 = zeros(1,L);
%% Model Parameters
% Initial weights
w_init.s_sp = 0.2;
w_init.d_sp = 0.2;
w_init.sp_l4 = 0.11;
w_init.s_l4 = 0.02;
w_init.d_l4 = 0.02;
w = w_init;
% Max value of Weights
maxw.s_l4 = 0.4;
maxw.d_l4 = 0.4;
maxw.sp_l4 = 0.11;
% Min Value of Weights
minw.s_l4 = 0.0001;
minw.d_l4 = 0.0001;
minw.sp_l4 = 0.0001;
% Short Term Plasticity
tau.th.re = 0.9;
tau.th.ei = 10;
tau.th.ir = tau_ir;
tau.sp.re = 0.9;
tau.sp.ei = 27;
tau.sp.ir = tau_ir;
x.s.e = 0;
x.s.r = 1;
x.s.i = 0;
x.d.e = 0;
x.d.r = 1;
x.d.i = 0;
x.sp.e = 0;
x.sp.r = 1;
x.sp.i = 0;
%
xarr.s.e = zeros(1,L);
xarr.s.r = zeros(1,L);
xarr.s.i = zeros(1,L);
xarr.d.e = zeros(1,L);
xarr.d.r = zeros(1,L);
xarr.d.i = zeros(1,L);
xarr.sp.e = zeros(1,L);
xarr.sp.r = zeros(1,L);
xarr.sp.i = zeros(1,L);
%
tau.syn = 10;
beta = 5;
tau.ref = 2;
%%
g.s = zeros(1,L);
g.d = zeros(1,L);
g.sp = zeros(1,L);
V.sp = zeros(1,L);
V.l4 = zeros(1,L);
%%
for t = 2:L
    if(Spikes.s(t-1))
        g.s(t:end) = g.s(t-1) + exp(-((t:L) - t)/tau.syn);
    elseif(Spikes.d(t-1))
        g.d(t:end) = g.d(t-1) + exp(-((t:L) - t)/tau.syn);
    elseif(Spikes.sp(t-1))
        g.sp(t:end) = g.sp(t-1) + exp(-((t:L) - t)/tau.syn);
    end
    V.sp(t) = V.sp(t) + g.s(t-1)*w.s_sp*x.s.e + g.d(t-1)*w.d_sp*x.d.e;
    V.l4(t) = V.l4(t) + g.s(t-1)*w.s_l4*x.s.e + g.d(t-1)*w.d_l4*x.d.e + g.sp(t-1)*w.sp_l4*x.sp.e;
    V.sp(t) = V.sp(t) - 0.1*abs(V.sp(t));
    V.l4(t) = V.l4(t) - 0.1*abs(V.l4(t));
    
    if(V.sp(t) > 0.05)
        Spikes.sp(t) = 1;
        V.sp(t+1:t+20) = V.sp(t) - beta*exp(-((t+1:t+20) - (t+1))/tau.ref);
    end
    
    if(V.l4(t) > 0.05)
        Spikes.l4(t) = 1;
        V.l4(t+1:t+20) = V.l4(t) - beta*exp(-((t+1:t+20) - (t+1))/tau.ref);
    end
    xarr.s.e(t-1) = x.s.e;
    xarr.s.r(t-1) = x.s.r;
    xarr.s.i(t-1) = x.s.i;
    xarr.d.e(t-1) = x.d.e;
    xarr.d.r(t-1) = x.d.r;
    xarr.d.i(t-1) = x.d.i;
    xarr.sp.e(t-1) = x.sp.e;
    xarr.sp.r(t-1) = x.sp.r;
    xarr.sp.i(t-1) = x.sp.i;
    x = ShortTermPlasticity(x,tau,Spikes.s(t),Spikes.d(t),Spikes.sp(t));
end
end
function [Input,L] = GenSpike()
p.ns = 0.5/1000*ones(1,250);
p.s = 10/1000*ones(1,50);
p.s1 = 2.5/1000*ones(1,50);
lamda.S = [p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s,p.ns...
    , p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s1,p.ns, p.s,p.ns, p.s,p.ns...
    , p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s,p.ns];
lamda.D = [p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s1,p.ns...
    , p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s,p.ns, p.s1,p.ns, p.s1,p.ns...
    , p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s1,p.ns];
L = length(lamda.S);
Input.s = zeros(1,L);
Input.d = zeros(1,length(lamda.D));
for j = 1:L
    Input.s(j) = binornd(1,lamda.S(j));
    Input.d(j) = binornd(1,lamda.D(j));
end
% subplot(2,1,1);
% plot(Input.s);
% title('S');
% xlabel('t(ms)');
% ylabel('Spikes');
% ylim([0 1]);
% subplot(2,1,2);
% plot(Input.d);
% title('D');
% xlabel('t(ms)');
% ylabel('Spikes');
% ylim([0 1]);
end
function y = ShortTermPlasticity(x,tau,s1,s2,s3)
y.s.r = x.s.r - s1*x.s.r/tau.th.re + x.s.i/tau.th.ir;
y.s.e = x.s.e + s1*x.s.r/tau.th.re - x.s.e/tau.th.ei;
y.s.i = x.s.i + x.s.e/tau.th.ei - x.s.i/tau.th.ir;

y.d.r = x.d.r - s2*x.d.r/tau.th.re + x.d.i/tau.th.ir;
y.d.e = x.d.e + s2*x.d.r/tau.th.re - x.d.e/tau.th.ei;
y.d.i = x.d.i + x.d.e/tau.th.ei - x.d.i/tau.th.ir;

y.sp.r = x.sp.r - s3*x.sp.r/tau.sp.re + x.sp.i/tau.sp.ir;
y.sp.e = x.sp.e + s3*x.sp.r/tau.sp.re - x.sp.e/tau.sp.ei;
y.sp.i = x.sp.i + x.sp.e/tau.sp.ei - x.sp.i/tau.sp.ir;
end
function plotspikes(Spikes)
figure;
subplot(4,1,1);
plot(Spikes.s);
xlabel('t(ms)');
ylabel('S');
ylim([0 1]);
subplot(4,1,2);
plot(Spikes.d);
xlabel('t(ms)');
ylabel('D');
ylim([0 1]);
subplot(4,1,3);
plot(Spikes.sp);
xlabel('t(ms)');
ylabel('SP');
ylim([0 1]);
subplot(4,1,4);
plot(Spikes.l4);
xlabel('t(ms)');
ylabel('L4');
ylim([0 1]);
end