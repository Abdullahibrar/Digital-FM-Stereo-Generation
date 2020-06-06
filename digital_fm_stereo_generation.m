clc;
clear;
close all;

noCh = 1; % number of channels, not important at the moment, just use 1
Fs = 400e3; % sampling frequency
bw_rx = 200e3; % receiver bandwidth (not being used atm)
i_t = 1; % integration time is always 1 second during FM signal generation
[data, Fs_data] = audioread('sound.mp3'); % music data
t = linspace(0, i_t, Fs*i_t); % time vector
st = zeros(noCh, Fs*i_t);
for k=1:noCh
    l = transpose(data(1e5+1:1e5+Fs_data, 1)); % left channel data
    r = transpose(data(1e5+1:1e5+Fs_data, 2)); % righ channel data

    l = resample(l, Fs*i_t, Fs_data*i_t); % interpolate to Fs for shifting in frequency domain
    r = resample(r, Fs*i_t, Fs_data*i_t); 
    l = l/max(l); % normalize
    r = r/max(r);

    figure
    subplot(2,1,1)

    % message signal
    mt = 0.5*(l+r) + 0.5*(l-r).*cos(2*pi*2*19000*t) + 0.5*cos(2*pi*19000*t);
    mt = mt/max(mt);
    temp = abs(fft(mt))/max(abs(fft(mt)));
    % plot the message signal
    plot(20*log10(fftshift(temp)))
    title('Message Signal');
    % integration and multiplication with 2*pi*75000
    mt = 2*pi*75000*cumsum(mt);
    st(noCh, :) = cos(mt) + 1i*sin(mt); % complex envelop stereo FM signal
end
subplot(2,1,2)
% plot the complex envelope signal
plot(20*log10(fftshift(abs(fft(st)))))
title('Complex Envelope Signal');