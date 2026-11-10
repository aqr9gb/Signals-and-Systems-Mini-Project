% =====================================================================
% Signals and Systems (3302-EE-3) - Mini Project - Project Set 10
% Tool    : GNU Octave
%
% Parts:  a) v(t)   b) w(t)   c) q(t)   d) r(t) = v(t) * w(t)
%         e) even and odd parts of r(t)   f) voice power spectrum
% Extra:  sampling-rate check for the convolution
%
% All figures are saved automatically as PNG files in the same folder.
% All printed numbers are also saved in the file results_log.txt.
% =====================================================================
clc; clear; close all;
more off;

% Save everything printed in the Command Window to a text file
if exist("results_log.txt", "file")
  delete("results_log.txt");
end
diary("results_log.txt");

% ---------------------------- Settings -------------------------------
fs = 44100;              % sampling rate (Hz)
dt = 1/fs;               % time between two samples (s)
T  = 8;                  % duration of the signals (s)
f  = 220;                % reference frequency (Hz)
t  = 0:dt:T;             % time axis from 0 to 8 s
N  = length(t);          % number of samples

zoom_end = 0.05;         % zoom window: first 0.05 s (11 cycles of 220 Hz)
fig_size = [0 0 9 6.5];  % size of the saved figures (inches)

% Unit step u(t): it is 1 for every t >= 0, so it is 1 on the whole axis
u = double(t >= 0);

printf("GNU Octave version %s\n", version);
printf("fs = %d Hz, dt = %.4e s, N = %d samples\n\n", fs, dt, N);

% ------------------------- Part a: v(t) ------------------------------
% Cosine of 220 Hz multiplied by a decaying exponential
v = cos(2*pi*f*t) .* exp(-1.3*t) .* u;

figure;
set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
         "paperposition", fig_size);
subplot(2,1,1);
plot(t, v); grid on;
xlabel("Time (s)"); ylabel("v(t)");
title("v(t) = cos(2\\pi 220 t) e^{-1.3t} u(t)  -  full range (0 to 8 s)");
subplot(2,1,2);
plot(t, v, "linewidth", 1.2); grid on; hold on;
plot(t,  exp(-1.3*t), "r--");            % upper envelope
plot(t, -exp(-1.3*t), "r--");            % lower envelope
xlim([0 zoom_end]);
xlabel("Time (s)"); ylabel("v(t)");
title("v(t)  -  zoom on the first 0.05 s  (red dashed lines: envelope \\pm e^{-1.3t})");
print("-dpng", "-r150", "fig1_v.png");

printf("Part a: max v = %.4f, min v = %.4f, v(8) = %.3e\n", ...
       max(v), min(v), v(end));

% ------------------------- Part b: w(t) ------------------------------
% Octave definition:  sinc(x) = sin(pi*x) / (pi*x),  and sinc(0) = 1
w = sinc(3*t) .* u;

figure;
set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
         "paperposition", fig_size);
subplot(2,1,1);
plot(t, w, "linewidth", 1.2); grid on;
xlabel("Time (s)"); ylabel("w(t)");
title("w(t) = sinc(3t) u(t)  -  full range (0 to 8 s)");
subplot(2,1,2);
plot(t, w, "linewidth", 1.2); grid on;
xlim([0 2]);
xlabel("Time (s)"); ylabel("w(t)");
title("w(t)  -  zoom on the first 2 s");
print("-dpng", "-r150", "fig2_w.png");

printf("Part b: max w = %.4f, min w = %.4f (at t = %.4f s)\n", ...
       max(w), min(w), t(w == min(w)));

% ------------------------- Part c: q(t) ------------------------------
q = 5 + v - w;

figure;
set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
         "paperposition", fig_size);
subplot(2,1,1);
plot(t, q); grid on;
xlabel("Time (s)"); ylabel("q(t)");
title("q(t) = 5 + v(t) - w(t)  -  full range (0 to 8 s)");
subplot(2,1,2);
plot(t, q, "linewidth", 1.2); grid on;
xlim([0 zoom_end]);
xlabel("Time (s)"); ylabel("q(t)");
title("q(t)  -  zoom on the first 0.05 s");
print("-dpng", "-r150", "fig3_q.png");

printf("Part c: q(0) = %.4f, max q = %.4f, min q = %.4f, q(8) = %.4f\n\n", ...
       q(1), max(q), min(q), q(end));

% --------------------- Part d: r(t) = v(t) * w(t) --------------------
% fftconv gives the discrete convolution (a sum). Multiplying by dt turns
% the sum into an approximation of the continuous convolution integral.
r_all = fftconv(v, w) * dt;      % length 2N-1, covers 0 to 16 s

% Only 0 to 8 s is valid: for t > 8 s the integral needs values of v and w
% after 8 s, but the signals were cut at 8 s.
r = r_all(1:N);

figure;
set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
         "paperposition", fig_size);
subplot(2,1,1);
plot(t, r); grid on;
xlabel("Time (s)"); ylabel("r(t)");
title("r(t) = v(t) * w(t)  -  valid range (0 to 8 s)");
subplot(2,1,2);
plot(t, r, "linewidth", 1.2); grid on;
xlim([0 zoom_end]);
xlabel("Time (s)"); ylabel("r(t)");
title("r(t)  -  zoom on the first 0.05 s");
print("-dpng", "-r150", "fig4_r.png");

% Measurements on r(t)
[r_peak, i_peak] = max(abs(r));
R      = abs(fft(r));
f_axis = (0:N-1) * fs / N;
[~, i_f] = max(R(1:floor(N/2)));
printf("Part d: peak |r| = %.4e at t = %.4f s\n", r_peak, t(i_peak));
printf("        main frequency of r = %.3f Hz\n", f_axis(i_f));
printf("        1/(2*pi*f) = %.4e (expected amplitude)\n", 1/(2*pi*f));
[~, i1] = min(abs(t - 1));
[~, i2] = min(abs(t - 2));
pk1 = max(abs(r(i1 - 200 : i1 + 200)));     % local peak of r near t = 1 s
pk2 = max(abs(r(i2 - 200 : i2 + 200)));     % local peak of r near t = 2 s
printf("        local peak near 1 s = %.4e, near 2 s = %.4e\n", pk1, pk2);
printf("        ratio = %.4f, exp(-1.3) = %.4f\n\n", pk2/pk1, exp(-1.3));

% ---------------- Part e: even and odd parts of r(t) -----------------
% Build a symmetric time axis from -8 to 8 s. r(t) is zero for t < 0.
t_sym = (-(N-1):(N-1)) * dt;
r_sym = [zeros(1, N-1), r];          % r(t) on the symmetric axis
r_rev = fliplr(r_sym);               % r(-t)

r_even = (r_sym + r_rev) / 2;        % even part
r_odd  = (r_sym - r_rev) / 2;        % odd part

figure;
set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
         "paperposition", fig_size);
subplot(2,1,1);
plot(t_sym, r_even); grid on;
xlim([-T T]); set(gca, "xtick", -T:2:T);
xlabel("Time (s)"); ylabel("r_e(t)");
title("Even part of r(t):  r_e(t) = [r(t) + r(-t)] / 2");
subplot(2,1,2);
plot(t_sym, r_odd); grid on;
xlim([-T T]); set(gca, "xtick", -T:2:T);
xlabel("Time (s)"); ylabel("r_o(t)");
title("Odd part of r(t):  r_o(t) = [r(t) - r(-t)] / 2");
print("-dpng", "-r150", "fig5_even_odd.png");

figure;
set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
         "paperposition", fig_size);
subplot(2,1,1);
plot(t_sym, r_even, "linewidth", 1.2); grid on;
xlim([-zoom_end zoom_end]);
xlabel("Time (s)"); ylabel("r_e(t)");
title("Even part  -  zoom around t = 0");
subplot(2,1,2);
plot(t_sym, r_odd, "linewidth", 1.2); grid on;
xlim([-zoom_end zoom_end]);
xlabel("Time (s)"); ylabel("r_o(t)");
title("Odd part  -  zoom around t = 0");
print("-dpng", "-r150", "fig6_even_odd_zoom.png");

% Numerical checks
pos = t_sym > 0;                                   % samples with t > 0
err_sum  = max(abs(r_even + r_odd - r_sym));       % even + odd = r ?
err_even = max(abs(r_even - fliplr(r_even)));      % r_e(t) =  r_e(-t) ?
err_odd  = max(abs(r_odd  + fliplr(r_odd)));       % r_o(t) = -r_o(-t) ?
err_half = max(abs(r_even(pos) - r_sym(pos)/2));   % r_e = r/2 for t > 0 ?
printf("Part e: max |r_e + r_o - r|       = %.3e\n", err_sum);
printf("        max |r_e(t) - r_e(-t)|    = %.3e\n", err_even);
printf("        max |r_o(t) + r_o(-t)|    = %.3e\n", err_odd);
printf("        max |r_e - r/2| for t > 0 = %.3e\n\n", err_half);

% ---------------- Part f: voice recording and FFT --------------------
voice_file = "voice.wav";

if exist(voice_file, "file")
  [x, fs_a] = audioread(voice_file);
  x  = mean(x, 2)';                 % one channel (mono), as a row
  Na = 8 * fs_a;                    % number of samples in exactly 8 s
  if length(x) >= Na
    x = x(1:Na);                    % keep the first 8 s
  else
    x = [x, zeros(1, Na - length(x))];   % fill with zeros up to 8 s
  end
  ta = (0:Na-1) / fs_a;

  figure;
  set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
           "paperposition", [0 0 9 4.5]);
  plot(ta, x); grid on;
  xlabel("Time (s)"); ylabel("Amplitude");
  title("Recorded voice signal (8 s)");
  print("-dpng", "-r150", "fig7_voice_time.png");

  % Power spectrum with the FFT
  X  = fft(x);
  P2 = abs(X).^2 / Na^2;            % two-sided power spectrum
  P1 = P2(1 : Na/2 + 1);            % keep 0 to fs/2
  P1(2:end-1) = 2 * P1(2:end-1);    % one-sided: double the inner bins
  fa = (0 : Na/2) * fs_a / Na;      % frequency axis (Hz)
  P_dB = 10 * log10(P1 + eps);      % power in dB

  % Smoothed spectrum: moving average over 41 bins (about 5 Hz).
  % Speech pitch changes a little with time, so smoothing makes the
  % peaks easier to read.
  M    = 41;
  P_s  = conv(P1, ones(1, M) / M, "same");
  Ps_dB = 10 * log10(P_s + eps);

  % Fundamental frequency: strongest smoothed peak between 70 and 300 Hz
  band = (fa >= 70) & (fa <= 300);           % usual range of voice pitch
  Pb   = P_s;  Pb(~band) = 0;
  [~, i0] = max(Pb);
  f0 = fa(i0);

  % Harmonics: strongest smoothed peak near each multiple of f0
  n_h = 6;
  f_h = zeros(1, n_h);
  L_h = zeros(1, n_h);
  for k = 1:n_h
    near = (fa >= k*f0 - f0/3) & (fa <= k*f0 + f0/3);
    Pk   = P_s;  Pk(~near) = 0;
    [~, ik] = max(Pk);
    f_h(k) = fa(ik);
    L_h(k) = Ps_dB(ik);
  end

  figure;
  set(gcf, "paperunits", "inches", "paperpositionmode", "manual", ...
           "paperposition", fig_size);
  subplot(2,1,1);
  plot(fa, P_dB); grid on;
  xlim([0 fs_a/2]);
  xlabel("Frequency (Hz)"); ylabel("Power (dB)");
  title("One-sided power spectrum of the voice (0 to fs/2)");
  subplot(2,1,2);
  plot(fa, P_dB, "color", [0.6 0.75 0.9]); grid on; hold on;
  plot(fa, Ps_dB, "b", "linewidth", 1.5);
  plot(f_h, L_h, "ro", "markersize", 8, "linewidth", 1.5);
  xlim([0 1000]); ylim([-110 -30]);
  xlabel("Frequency (Hz)"); ylabel("Power (dB)");
  title(["Zoom on 0 to 1000 Hz  (dark line: smoothed spectrum, ", ...
         "red circles: f_0 and its harmonics)"]);
  print("-dpng", "-r150", "fig8_voice_spectrum.png");

  % Measurements on the voice
  printf("Part f: fs = %d Hz, samples = %d, resolution = %.3f Hz\n", ...
         fs_a, Na, fs_a/Na);
  printf("        peak amplitude = %.4f\n", max(abs(x)));
  printf("        mean power (time)   = %.4e\n", mean(x.^2));
  printf("        sum of the spectrum = %.4e\n", sum(P1));
  printf("        fundamental frequency f0 = %.1f Hz\n", f0);
  for k = 1:n_h
    printf("        peak %d: %6.1f Hz = %.2f x f0, level = %6.1f dB\n", ...
           k, f_h(k), f_h(k)/f0, L_h(k));
  end
  P_tot = sum(P1(2:end));                    % total power without DC
  edges = [300 1000 4000 8000];
  for k = 1:length(edges)
    part = sum(P1(fa > 0 & fa <= edges(k))) / P_tot * 100;
    printf("        power below %4d Hz = %.2f %%\n", edges(k), part);
  end
else
  printf("Part f skipped: file %s was not found in this folder.\n", voice_file);
end

% ------------- Extra: sampling-rate check for the convolution --------
% r(t) is computed again with the same method (fftconv multiplied by dt)
% at three sampling rates. Each result is compared with a reference
% computed at a very high sampling rate. The check uses the first 1 s,
% where r(t) is largest. This is allowed because r(t) at time t depends
% only on v and w between 0 and t.
T_test = 1;                   % test duration (s)
fs_ref = 3528000;             % reference rate: a multiple of all test rates
fs_all = [2000 8000 44100];   % sampling rates to test (Hz)

t_ref = 0 : 1/fs_ref : T_test;
v_ref = cos(2*pi*f*t_ref) .* exp(-1.3*t_ref) .* double(t_ref >= 0);
w_ref = sinc(3*t_ref) .* double(t_ref >= 0);
r_ref = fftconv(v_ref, w_ref) / fs_ref;
r_ref = r_ref(1 : length(t_ref));

printf("\nSampling-rate check (first %d s):\n", T_test);
printf("  reference: fs = %d Hz, peak |r| = %.4e\n", fs_ref, max(abs(r_ref)));

for k = 1:length(fs_all)
  fs_k = fs_all(k);
  t_k  = 0 : 1/fs_k : T_test;
  v_k  = cos(2*pi*f*t_k) .* exp(-1.3*t_k) .* double(t_k >= 0);
  w_k  = sinc(3*t_k) .* double(t_k >= 0);
  r_k  = fftconv(v_k, w_k) / fs_k;
  r_k  = r_k(1 : length(t_k));

  step  = fs_ref / fs_k;               % reference samples per test sample
  r_cmp = r_ref(1 : step : end);       % reference at the same time points
  err_k = max(abs(r_k - r_cmp)) / max(abs(r_ref)) * 100;
  printf("  fs = %6d Hz : peak |r| = %.4e , max error = %6.2f %% of the peak\n", ...
         fs_k, max(abs(r_k)), err_k);
end

printf("\nDone. All figures were saved as PNG files.\n");
diary off;
