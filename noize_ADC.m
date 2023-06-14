% digital filter 
% 32bit signal 

% noize 4ksapmle 16uV
N = 1e4;
noize_level = 8e-6;        % level noize
f_dateRate = 4e3;           % f_dateRate




t = 1/f_dateRate:1/f_dateRate:N/f_dateRate;

Gause_Voltage = (noize_level/2).*wgn(N,1,0); 


Gause_Voltage = Gause_Voltage';% + sin(2*pi*100.*t);

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

%% graph median Filter impulse pomeh
figure();
vect_signal_median = medfilt1(Gause_Voltage,9);

Wind = 5;

for i = 1 : N
    if (i > floor(Wind/2) && i <= N-floor(Wind/2))
        sig_sort = sort_buble(Gause_Voltage(i-floor(Wind/2):(i+floor(Wind/2))), Wind);
        vect_signal_med(i) = sig_sort(round(Wind/2));
        summ = 0;
        buff = Gause_Voltage(i-floor(Wind/2):(i+floor(Wind/2)));
        for j = 1 : Wind
            summ = summ + buff(j);
        end;
        Moving_Aver(i) = summ/Wind;
    else
        vect_signal_med(i) = 0;
        Moving_Aver(i) = 0;
    end;
end;

plot(t, Gause_Voltage, t, vect_signal_median, t, vect_signal_med, t, Moving_Aver); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize median filter model');

%% graph spectr
figure();
S = fft(vect_signal_med);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, abs(S)); grid;
title("Spectr noize median filter fd=4kHz");

%% Gausse Filter
W = Num;

voltage_Gause = filter(W, 1, Gause_Voltage);
figure();

plot(t, Gause_Voltage, t, voltage_Gause); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize Gausse Filter');

%% graph spectr
figure();
S = fft(voltage_Gause);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, abs(S)); grid;
title("Spectr noize Window filter fd=4kHz");

%% filter alpha-betta Kalman 
figure();
dt = 0.5;
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

plot(t, Gause_Voltage, t, voltage_noze_Kalman); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize Kalman alpha-betta Filter');

%% graph spectr
figure();
S = fft(voltage_noze_Kalman);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, abs(S)); grid;
title("Spectr noize Kalman Alpha-Betta Filter fd=4kHz");

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