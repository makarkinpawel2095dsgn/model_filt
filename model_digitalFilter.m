% digital filter 
% 32bit signal 

% noize 4ksapmle 16uV
N = 1e4;
noize_level = 50e-6;        % level noize
f_dateRate = 4e3;           % f_dateRate




t = 1/f_dateRate:1/f_dateRate:N/f_dateRate;

Gause_Voltage = (noize_level/4).*wgn(N,1,0); 


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
plot(df, abs(S)); grid;
title("Spectr noize fd=4kHz");

%% graph median Filter impulse pomeh
figure();
vect_signal_median = medfilt1(Gause_Voltage,17);
plot(t, vect_signal_median); grid;
ylabel('Voltage, uV');
xlabel('Time, s');
title('Signal Noize median filter model');

%% graph spectr
figure();
S = fft(vect_signal_median);
df = f_dateRate/N : f_dateRate/N : f_dateRate;
plot(df, abs(S)); grid;
title("Spectr noize fd=4kHz");


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

