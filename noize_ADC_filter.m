% digital filter 
% 32bit signal 
clear;
clc;
% noize 4ksapmle 16uV
N = 1e4;
noize_level = 6e-6;        % level noize
f_dateRate = 4e3;           % f_dateRate

t = 1/f_dateRate:1/f_dateRate:N/f_dateRate;

Gause_Voltage = (noize_level/2).*wgn(N,1,0); 

Gause_Voltage = Gause_Voltage';% + sin(2*pi*10.*t);

%% graph signal noize
figure();
plot(t, Gause_Voltage); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize model');


%% graph spectr
figure();
S = fft(Gause_Voltage);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, 20.*log10(abs(S))); grid;
title("Spectr noize fd=4kHz");


%% my Filter
% Coeff = FIR_FILTER_LPF(500, 4000, 60, 16);
Coeff = FIR_FILTER_HPF(1000, 4000, 60, 64);
out_FIR_HPF = FIR_FILTER(Gause_Voltage, Coeff);
out_Proc = Gause_Voltage - out_FIR_HPF;
Coeff = FIR_FILTER_LPF(800, 4000, 60, 32);
out_FIR_LPF = FIR_FILTER(out_Proc, Coeff);
figure();
plot(t, Gause_Voltage, t, out_FIR_LPF); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize filter model');
figure();
S = fft(out_FIR_LPF);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, 20.*log10(abs(S))); grid;
title("Spectr noize my filter fd=4kHz");
%% graph median Filter impulse pomeh
figure();
vect_signal_median = medfilt1(Gause_Voltage,9);

Wind = 101;

for i = 1 : N
    if (i > floor(Wind/2) && i <= N-floor(Wind/2))
        sig_sort = sort_buble(Gause_Voltage(i-floor(Wind/2):(i+floor(Wind/2))), Wind);
        
        std_array = stdArrow(Gause_Voltage(i-floor(Wind/2):(i+floor(Wind/2))), Wind);
        corr_vect = auto_corr(Gause_Voltage(i-floor(Wind/2):(i+floor(Wind/2))), Wind);
        
        vect_signal_med(i) = sig_sort(round(Wind/2));
        summ = 0;
        buff = Gause_Voltage(i-floor(Wind/2):(i+floor(Wind/2)));
        for j = 1 : Wind
            summ = summ + buff(j);
        end;
        Moving_Aver(i) = summ/Wind;
        stdFilter(i) = Gause_Voltage(i) - std_array;
    else
        vect_signal_med(i) = 0;
        Moving_Aver(i) = 0;
        stdFilter(i) = 0;
    end;
end;

out_FIR = FIR_FILTER(Gause_Voltage, [1 1 1 1 1 1 1 1]);

plot(t, Gause_Voltage, t, vect_signal_median, t, vect_signal_med, t, out_FIR); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize median filter model');

%% graph spectr
figure();
S = fft(out_FIR);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, abs(S)); grid;
title("Spectr noize median filter fd=4kHz");

% %% Gausse Filter
% W = Num;
% 
% voltage_Gause = filter(W, 1, Gause_Voltage);
% figure();
% 
% plot(t, Gause_Voltage, t, voltage_Gause); grid;
% ylabel('Voltage, uV');
% xlabel('Time, s');
% title('Signal Noize Gausse Filter');
% 
% %% graph spectr
% figure();
% S = fft(voltage_Gause);
% df = f_dateRate/N : f_dateRate/N : f_dateRate;
% plot(df, abs(S)); grid;
% title("Spectr noize Window filter fd=4kHz");

%% filter alpha-betta Kalman 
figure();
dt = 0.4;
alpha = 0.1;
betta = 0.005;

voltage_noze_to = 0;
v_to = 0;

for i = 1 : N

    voltage_noze = voltage_noze_to + (v_to*dt);
    voltage_k = v_to;

    rk = Gause_Voltage(i) - voltage_noze;

    voltage_noze = voltage_noze + alpha*rk;
    voltage_k = voltage_k + (betta*rk)/dt;

    voltage_noze_to = voltage_noze;
    v_to = voltage_k;

    voltage_noze_Kalman(i) = voltage_noze_to;
end;

dr = 0.05; dvx = 1/4e3; T = 0.1; 
qx = 1e-5; qvx = 1e-3;

H = [1 0]; F = [1 T; 0 1]; Q = [qx 0; 0 qvx]; R = dr;

xe = [Gause_Voltage(1); 0];  xf = [Gause_Voltage(1); 0];

De = [dr 0; 0 dvx];
D   = [dr 0; 0 dvx];

aXest = zeros(N,1);

for i=1:N
    xe = F * xf;
    De = F * D * F' + Q;
    
    S  = H*De*H' + R;
    
    K = De * H' / ( H*De*H' + R);
    
    xf = xe + K * (Gause_Voltage(i) - H * xe);
    
    D = ([1 0; 0 1] - K * H) * De;
    
    Kalman_lineir(i) = xf(1);
    
end

subplot(2,1,1); plot(t, Gause_Voltage, t, voltage_noze_Kalman); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize Kalman alpha-betta Filter');
subplot(2,1,2); plot(t, Gause_Voltage, t, Kalman_lineir); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize Kalman Filter');

%% graph spectr
figure();
S = fft(voltage_noze_Kalman);
S0 = fft(Kalman_lineir);

df = f_dateRate/N : f_dateRate/N : f_dateRate;
subplot(2,1,1); plot(df, abs(S)); grid;
title("Spectr noize Kalman Alpha-Betta Filter fd=4kHz");
subplot(2,1,2); plot(df, abs(S0)); grid;

%title("Spectr noize Kalman Alpha-Betta Filter fd=4kHz");

% %% Autocorrelation algorithm 
% 
% M = 500;
% figure();
% for i = 1 : N
%     if (i > floor(M/2) && i <= N-floor(M/2))
%         
%         buff = Gause_Voltage(i-floor(M/2):(i+floor(M/2)));
%         
%         summ= zeros(1,2*M);
%         for j = 1 : M
%             
%             for ij = 1 : M
%                 summ(j+ij-1) = summ(j+ij-1) + buff(j)*buff(ij);
%             end;
%            
%         end;
%         Autocorr = summ;
% 
%     end;
% end;
% plot(Autocorr);

function sig_sort = sort_buble(buff_buble, Window)
% signal - signal
% Window - window sort buble

for i = 1 : Window
    for j = 1 : Window-i
        %
        if ( buff_buble(j) > buff_buble(j+1) )
            var_buff = buff_buble(j);
            buff_buble(j) = buff_buble(j+1);
            buff_buble(j+1) = var_buff;
        end;
    end;
end;
sig_sort = buff_buble;
end


function std_array = stdArrow(buff_buble, Window)
% signal - signal
% Window - window sort buble
mean = 0;
for i = 1 : Window
    mean = mean + buff_buble(i);
end;

mean = mean/Window;

stds = 0;
for i = 1 : Window
    stds = stds + (mean - buff_buble(i))^2;
end;
std_array = stds/Window;
end


function corr_vect = auto_corr(buff_buble, Window)
% signal - signal
% Window - window sort buble
summ = zeros(1, Window);
for i = 1 : Window
    summ(i) = 0;
    for j = 1 : Window
        if (i-j>0) 
            summ(i) = summ(i) + buff_buble(j)*buff_buble(i-j);
        end;
    end;
end;
corr_vect = summ/Window;
end

function out_FIR = FIR_FILTER(in_FIR, Windw)
% signal - signal
% Window - window sort buble
summ = zeros(1, length(in_FIR));
for i = 1 : length(in_FIR)
    summ(i) = 0;
    for j = 1 : length(Windw)
        if (i-j>0) 
            summ(i) = summ(i) + Windw(j)*in_FIR(i-j);
        end;
    end;
end
out_FIR = summ/length(Windw);
end

function Coeff = FIR_FILTER_LPF(fs, fd, G_dB, N_FIR)
%%
% fs - freqwency fs - 0.7dB
% fd - freqwency f sample 
% Ap - 1;
% An - 0;
N_Fs = fs/fd;
N_Fp = (fd-fs)/fd;
% ones vect
Wind_one = ones(round(N_Fs*N_FIR/2), 1);
% null vect
Wind_null = zeros(round(N_Fp*N_FIR/2), 1);

Amp_g = 10^(-G_dB/20);

for i = 1 : length(Wind_null)
    Wind_null(i) = Amp_g;
end;

Wind = [Wind_one; Wind_null];
Window = [Wind' rot90(rot90(Wind'))];
Coeff = fft(Window);
Coeff = real([Coeff(N_FIR/2+1:N_FIR), Coeff(1:N_FIR/2)]);

Coeff = Coeff;%./sum(Coeff);
stem(Coeff);

title('Coeff FIR Filter');
[h, w] = freqz(Coeff, 1, 512);
figure
plot(w/(pi),20*log10(abs(h)/max(abs(h)))); grid;
%% window FIR
Wind_comp_Hemming = 0.54 - 0.46.*cos(2*pi*(1:length(Coeff))./length(Coeff));
Coeff_Hemming = Coeff.*Wind_comp_Hemming;
figure();
stem(Coeff_Hemming); grid;
title('Coeff FIR LPF Window Filter');
[h, w] = freqz(Coeff_Hemming, 1, 512);
figure
plot(w/(pi),20*log10(abs(h)/max(abs(h)))); grid;
Coeff = Coeff_Hemming;
end



function Coeff = FIR_FILTER_HPF(fs, fd, G_dB, N_FIR)
%%
% fs - freqwency fs - 0.7dB
% fd - freqwency f sample 
% Ap - 1;
% An - 0;
N_Fs = fs/fd;
N_Fp = (fd-fs)/fd;
% ones vect
Wind_one = ones(round(N_Fp*N_FIR/2), 1);
% null vect
Wind_null = zeros(round(N_Fs*N_FIR/2), 1);

Amp_g = 10^(-G_dB/20);

for i = 1 : length(Wind_null)
    Wind_null(i) = Amp_g;
end;

Wind = [Wind_null; Wind_one];
Window = [Wind' rot90(rot90(Wind'))];
Coeff = fft(Window);
Coeff = real([Coeff(N_FIR/2+1:N_FIR), Coeff(1:N_FIR/2)]);

Coeff = Coeff;%./sum(Coeff);
stem(Coeff); grid;
title('Coeff FIR LPF Filter');
[h, w] = freqz(Coeff, 1, 512);
figure
plot(w/(pi),20*log10(abs(h)/max(abs(h)))); grid;

%% window FIR
Wind_comp_Hemming = 0.54 - 0.46.*cos(2*pi*(1:length(Coeff))./length(Coeff));
Coeff_Hemming = Coeff.*Wind_comp_Hemming;
figure();
stem(Coeff_Hemming); grid;
title('Coeff FIR LPF Window Filter');
[h, w] = freqz(Coeff_Hemming, 1, 512);
figure
plot(w/(pi),20*log10(abs(h)/max(abs(h)))); grid;

Coeff = Coeff_Hemming;
end