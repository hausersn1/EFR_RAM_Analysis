%% RAM Generator

fs = 48828; % Sampling rate
duration = 1.0; % Second
fm = 223; % Modulation frequency
fc = 4000; % Carrier frequency

t = 0:(1/fs):(duration - 1/fs);
carr = sin(2*pi*fc*t);

mod = (square(2*pi*fm*t, 25) + 1)/2;
signal = (mod .* carr)*0.1;

signal_rms = rms(signal);
y = scaleSound(rampsound(signal, fs, 0.01));

audiowrite('RAM_223.wav', y, fs); 
