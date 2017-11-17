clear;
clc;
% %%
% n = 700;
% lamda  = 40+rand(1,1000)*(50-40);
% p = lamda./1000;         %probabilty of spiking per 1 ms
% Spikes = zeros(n,1000);
% for i = 1:n
%     for j = 1:1000
%         if(rand() < p(j))
%             Spikes(i,j) = 1;
%         end
%     end
% end
% SpikeTime = cell(320,1);
% for i = 1:n
%     SpikeTime{i,1} = find(Spikes(i,:));
% end
% %%
% a = [10 20 40 80 160 320 700];
% mse = zeros(1,length(a));
% PSTH = zeros(length(a),1000);
% for i = 1:length(a)
%     PSTH(i,:) = 1000*mean(Spikes(1:a(i),:),1);
%     mse(i) = sqrt(immse(lamda,PSTH(i,:)));
% end
% plot(a,mse,'x-');
% ylabel('Mean Square Error');
% xlabel('# of Repetitions');
% figure;
%%
% clearvars p lamda;
p.ns = 0.5/1000*ones(1,250);
p.s = 10/1000*ones(1,50);
p.s1 = 2.5/1000*ones(1,50);
lamda.S = [p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s,p.ns...
    , p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s1,p.ns, p.s,p.ns, p.s,p.ns...
    , p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s,p.ns, p.s,p.ns];
lamda.D = [p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s1,p.ns...
    , p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s,p.ns, p.s1,p.ns, p.s1,p.ns...
    , p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s1,p.ns, p.s1,p.ns];
Input.s = zeros(1,length(lamda.S));
Input.d = zeros(1,length(lamda.D));
for j = 1:length(lamda.S)
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
Spikes = Input;
Spikes.sp = zeros(1,length(lamda.S));
Spikes.l4 = zeros(1,length(lamda.S));
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
tau.th.ir = 5000;
tau.sp.re = 0.9;
tau.sp.ei = 27;
tau.sp.ir = 5000;
x.s.e = 1;
x.s.r = 0;
x.s.i = 0;
x.d.e = 1;
x.d.r = 0;
x.d.i = 0;
x.sp.e = 1;
x.sp.r = 0;
x.sp.i = 0;
%
tau.syn = 10;
beta = 5;
tau.ref = 2;
%%
g.s = zeros(1,length(lamda.S));
g.d = zeros(1,length(lamda.S));
g.sp = zeros(1,length(lamda.S));
V.sp = zeros(1,length(lamda.S));
V.l4 = zeros(1,length(lamda.S));
%%
for t = 2:length(lamda.S)
    if(Spikes.s(t-1))
        g.s(t:end) = g.s(t-1) + exp(-((t:length(lamda.S)) - t)/tau.syn);
    elseif(Spikes.d(t-1))
        g.d(t:end) = g.d(t-1) + exp(-((t:length(lamda.S)) - t)/tau.syn);
    elseif(Spikes.sp(t-1))
        g.sp(t:end) = g.sp(t-1) + exp(-((t:length(lamda.S)) - t)/tau.syn);
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
    x = ShortTermPlasticity(x,tau,Spikes.s(t),Spikes.d(t),Spikes.sp(t));
end
figure(3);
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
% figure(4);
% subplot(2,1,1);
% plot(V.sp);
% xlabel('t(ms)');
% ylabel('V_SP');
% subplot(2,1,2);
% plot(V.l4);
% xlabel('t(ms)');
% ylabel('V_L4');
% 
% figure(5);
% plot(g.s);
% xlabel('t(ms)');
% ylabel('g.s');
function y = ShortTermPlasticity(x,tau,s1,s2,s3)
y.s.r = x.s.r - s1*x.s.r/tau.th.re + x.s.i/tau.th.ir;
y.s.e = x.s.e + s1*x.s.r/tau.th.re - x.s.e/tau.th.ei;
y.s.i = x.s.e/tau.th.ei - x.s.i/tau.th.ir;

y.d.r = x.d.r - s2*x.d.r/tau.th.re + x.d.i/tau.th.ir;
y.d.e = x.d.e + s2*x.d.r/tau.th.re - x.d.e/tau.th.ei;
y.d.i = x.d.e/tau.th.ei - x.d.i/tau.th.ir;

y.sp.r = x.sp.r - s3*x.sp.r/tau.sp.re + x.sp.i/tau.sp.ir;
y.sp.e = x.sp.e + s3*x.sp.r/tau.sp.re - x.sp.e/tau.sp.ei;
y.sp.i = x.sp.e/tau.sp.ei - x.sp.i/tau.sp.ir;
end