// ====================================================================
// CASE STUDY: Noise & Modulation Impact on Wireless Signals
// LANGUAGE: Scilab
// ====================================================================

clc;
clear;
rand("seed", 2025);

// ------------------ PARAMETERS ------------------
EbN0dB = 0:2:20;        // Eb/N0 values (dB)
Nbits = 1e5;            // Number of bits
BER_BPSK = zeros(1, length(EbN0dB));
BER_QPSK = zeros(1, length(EbN0dB));
BER_16QAM = zeros(1, length(EbN0dB));


// ------------------ FUNCTION ------------------
// Gray mapping 2 bits → PAM level (-3,-1,+1,+3)
function level = gray_map_2bit(b1, b2)
    if b1==0 & b2==0 then
        level = -3;
    elseif b1==0 & b2==1 then
        level = -1;
    elseif b1==1 & b2==1 then
        level = 1;
    else
        level = 3;
    end
endfunction

// ------------------ BPSK SIMULATION ------------------
disp("Simulating BPSK...");
for i = 1:length(EbN0dB)
    bits = round(rand(1, Nbits));
    s = 2*bits - 1;                        // 0→-1, 1→+1
    EbN0 = 10^(EbN0dB(i)/10);
    N0 = 1/EbN0;
    noise = sqrt(N0/2) * rand(1, Nbits, "normal");
    r = s + noise;
    bits_hat = (r >= 0);
    BER_BPSK(i) = sum(bits <> bits_hat) / Nbits;
    printf("BPSK Eb/N0=%2d dB -> BER=%1.4e\n", EbN0dB(i), BER_BPSK(i));
end

// ------------------ QPSK SIMULATION ------------------
disp("Simulating QPSK...");
for i = 1:length(EbN0dB)
    Nsym = Nbits / 2;
    bits = round(rand(1, 2*Nsym));

    sI = 2*bits(1:2:2*Nsym) - 1;
    sQ = 2*bits(2:2:2*Nsym) - 1;
    s = (sI + %i*sQ) / sqrt(2);

    EbN0 = 10^(EbN0dB(i)/10);
    N0 = 1/EbN0;
    noise = sqrt(N0/2)*(rand(1,Nsym,"normal") + %i*rand(1,Nsym,"normal"));
    r = s + noise;

    bits_hatI = (real(r) >= 0);
    bits_hatQ = (imag(r) >= 0);

    bits_hat = zeros(1, 2*Nsym);
    bits_hat(1:2:2*Nsym) = bits_hatI;
    bits_hat(2:2:2*Nsym) = bits_hatQ;

    BER_QPSK(i) = sum(bits <> bits_hat) / Nbits;
    printf("QPSK Eb/N0=%2d dB -> BER=%1.4e\n", EbN0dB(i), BER_QPSK(i));
end

// ------------------ 16-QAM SIMULATION ------------------
disp("Simulating 16-QAM...");
for i = 1:length(EbN0dB)
    Nsym = Nbits / 4;
    bits = round(rand(1, 4*Nsym));
    sI = zeros(1, Nsym);
    sQ = zeros(1, Nsym);

    for k = 1:Nsym
        sI(k) = gray_map_2bit(bits(4*k-3), bits(4*k-2));
        sQ(k) = gray_map_2bit(bits(4*k-1), bits(4*k));
    end

    s = (sI + %i*sQ) / sqrt(10); // normalization
    EbN0 = 10^(EbN0dB(i)/10);
    N0 = 1/EbN0;
    noise = sqrt(N0/2)*(rand(1,Nsym,"normal")+%i*rand(1,Nsym,"normal"));
    r = s + noise;

    levels = [-3, -1, 1, 3];
    gray = [0 0; 0 1; 1 1; 1 0];
    bits_hat = zeros(1, 4*Nsym);

    for k = 1:Nsym
        diffI = abs(real(r(k))*sqrt(10) - levels);
        diffQ = abs(imag(r(k))*sqrt(10) - levels);
        [valI, idxI] = min(diffI);
        [valQ, idxQ] = min(diffQ);
        bits_hat(4*k-3:4*k-2) = gray(idxI, :);
        bits_hat(4*k-1:4*k) = gray(idxQ, :);
    end

    BER_16QAM(i) = sum(bits <> bits_hat) / Nbits;
    printf("16-QAM Eb/N0=%2d dB -> BER=%1.4e\n", EbN0dB(i), BER_16QAM(i));
end

// ------------------ BER TABLE ------------------
disp("-----------------------------------------------------");
disp(" Eb/N0(dB)    BER_BPSK       BER_QPSK       BER_16QAM ");
disp("-----------------------------------------------------");
for i = 1:length(EbN0dB)
    mprintf("   %2d           %1.4e      %1.4e      %1.4e\n", EbN0dB(i), BER_BPSK(i), BER_QPSK(i), BER_16QAM(i));
end
disp("-----------------------------------------------------");

// ------------------ LIVE BER PLOT ------------------
// keep this exactly as your working code, using scf(1)
scf(1);
for i = 1:length(EbN0dB)
    clf(); // clears figure 1 only
    semilogy(EbN0dB(1:i), BER_BPSK(1:i), 'b-o');
    plot(EbN0dB(1:i), BER_QPSK(1:i), 'r-s');
    plot(EbN0dB(1:i), BER_16QAM(1:i), 'g-^');
    xlabel("Eb/N0 (dB)");
    ylabel("Bit Error Rate (BER)");
    title("Live BER vs Eb/N0 for BPSK, QPSK, 16-QAM");
    legend(["BPSK", "QPSK", "16-QAM"]);
    xgrid();
    sleep(300);
end

// ------------------ CONSTELLATION ANIMATION ------------------
disp("Displaying live QPSK constellation...");
scf(2); // constellation in figure 2 (will not clash with waveform figs)
for i = 1:length(EbN0dB)
    Nsym = 1000;
    bits = round(rand(1, 2*Nsym));
    sI = 2*bits(1:2:2*Nsym) - 1;
    sQ = 2*bits(2:2:2*Nsym) - 1;
    s = (sI + %i*sQ) / sqrt(2);
    EbN0 = 10^(EbN0dB(i)/10);
    N0 = 1/EbN0;
    noise = sqrt(N0/2)*(rand(1,Nsym,"normal")+%i*rand(1,Nsym,"normal"));
    r = s + noise;

    clf(); // clears figure 2 only
    subplot(1,2,1);
    plot(real(s), imag(s), 'bo');
    title("Original QPSK Constellation");
    xlabel("In-Phase"); ylabel("Quadrature");
    xgrid();

    subplot(1,2,2);
    plot(real(r), imag(r), 'r.');
    title(msprintf("Noisy QPSK (Eb/N0 = %d dB)", EbN0dB(i)));
    xlabel("In-Phase"); ylabel("Quadrature");
    xgrid();

    sleep(400);
end

// ------------------ INTERACTIVE MODE ------------------
disp("-----------------------------------------------------");
disp("Interactive mode: visualize any Eb/N0 value you choose!");
SNR = input("Enter Eb/N0 value (dB): ");

Nsym = 1000;
bits = round(rand(1, 2*Nsym));
sI = 2*bits(1:2:2*Nsym) - 1;
sQ = 2*bits(2:2:2*Nsym) - 1;
s = (sI + %i*sQ)/sqrt(2);
EbN0 = 10^(SNR/10);
N0 = 1/EbN0;
noise = sqrt(N0/2)*(rand(1,Nsym,"normal")+%i*rand(1,Nsym,"normal"));
r = s + noise;

scf(3); // interactive constellation in figure 3
clf();
plot(real(r), imag(r), 'r.');
title(msprintf("QPSK Constellation at Eb/N0 = %d dB", SNR));
xlabel("In-Phase"); ylabel("Quadrature");
xgrid();

// ------------------ WAVEFORM VISUALIZATION ------------------
disp("Displaying waveform impact for BPSK...");

// keep Nbits_wave small and use a unique figure number (4)
Nbits_wave = 100;   // number of bits to visualize
bits_wave = round(rand(1, Nbits_wave));
s_wave = 2*bits_wave - 1;     // BPSK symbols
EbN0dB_wave = 6;    // example SNR
EbN0 = 10^(EbN0dB_wave/10);
N0 = 1/EbN0;
noise_wave = sqrt(N0/2) * rand(1, Nbits_wave, "normal");
r_wave = s_wave + noise_wave;

// Create a "time axis" (just for plotting)
t_wave = 1:Nbits_wave;

scf(4); // DISCRETE-TIME waveform in figure 4
clf();
subplot(2,1,1);
plot(t_wave, s_wave, 'b-', 'LineWidth', 2);
title("Original Transmitted BPSK Signal (Discrete Time)");
xlabel("Bit Index");
ylabel("Amplitude");
xgrid();

subplot(2,1,2);
plot(t_wave, r_wave, 'r-', 'LineWidth', 2);
title(msprintf("Received Noisy BPSK Signal (Eb/N0 = %d dB)", EbN0dB_wave));
xlabel("Bit Index");
ylabel("Amplitude");
xgrid();
sleep(800); // give time to view

// ------------------ CONTINUOUS CARRIER BPSK WAVEFORM ------------------
disp("Displaying continuous carrier-modulated BPSK waveform...");

// Use a unique figure number (5) so previous figs remain
fc = 1000;       // carrier frequency (Hz)
fs = 10000;      // sampling frequency (Hz)
Tb = 0.001;      // bit duration (1 ms)

// ensure fs*Tb is integer
samples_per_bit = fs*Tb;
if floor(samples_per_bit) <> samples_per_bit then
    // choose nearest integer and warn
    samples_per_bit = floor(samples_per_bit);
    warning("fs*Tb was not integer — using floor(samples_per_bit) to avoid indexing errors.");
end

t_cont = 0:1/fs:Nbits_wave*Tb - 1/fs;

carrier = cos(2*%pi*fc*t_cont);
waveform = zeros(1, length(t_cont));

for k = 1:Nbits_wave
    idx = (k-1)*samples_per_bit + 1 : k*samples_per_bit;
    // ensure idx doesn't exceed waveform length
    idx(idx > length(waveform)) = [];
    waveform(idx) = (2*bits_wave(k)-1) * carrier(idx);
end

noise_cont = sqrt(N0/2) * rand(1, length(t_cont), "normal");
r_wave_cont = waveform + noise_cont;

scf(5); // ANALOG carrier in figure 5
clf();
subplot(2,1,1);
plot(t_cont(1:min(2000,length(t_cont))), waveform(1:min(2000,length(t_cont)))); // show first few ms
title("Clean BPSK Carrier Signal");
xlabel("Time (s)");
ylabel("Amplitude");
xgrid();

subplot(2,1,2);
plot(t_cont(1:min(2000,length(t_cont))), r_wave_cont(1:min(2000,length(t_cont))));
title(msprintf("Noisy BPSK Carrier Signal (Eb/N0 = %d dB)", EbN0dB_wave));
xlabel("Time (s)");
ylabel("Amplitude");
xgrid();
sleep(1200); // give more time to inspect carrier plots

disp("-----------------------------------------------------");
disp("Simulation completed successfully!");
