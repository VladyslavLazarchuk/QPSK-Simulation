% Vladyslav Lazarchuk
% April 28 2025

% ECE644 Wireless Communications: From Fundamentals to 5G
% Dr. Joerg Kliewer
% Midterm Project: QPSK Matlab Simulation

clear
rng(9900) %Using the last 4 digits of my ID as a seed

    % Generating a sequence of random bits
bits =  randi([0,1], 1, 1e6);

    % Since QPSK takes 2 bits per symbol - we have to group the bits
symbols = reshape(bits, 2, []).';
map = grayMap(symbols); % Actual mapping

    % Generating complex-values AWGN samples:
samples = length(bits)/2; %number of samples
nR = randn(samples, 1); %Real part of AWGN
nI = randn(samples, 1); %Imaginary part of AWGN
noise = nR + nI*1i; %Combined AWGN

    %Initializing
ebN0_range = 0:8;
    %Initializing arrays to store data for cumulative plotting
BER = NaN(size(ebN0_range)); 
TBER = NaN(size(ebN0_range));

    %Plotting
figure;
    %Defining axis
xlim([0 8]);
ylim([10^-4 .1]);
yscale('log');
hold on;
    %Empty plots that will help us connect discrete values
empiricalPlot = semilogy(NaN, NaN, 'b-.', 'LineWidth', 1);
theoreticalPlot = semilogy(NaN, NaN, 'r-.', 'LineWidth', 1);

for ebN0_dB = ebN0_range %Loop that will go through all integer dB values
    nas = getNoiseAndSignal(ebN0_dB, noise, map); %Getting a signal mixed with AWGN
    recovSignal = remap(nas); %Gray mapping received signal
    recovBits = demodulate(recovSignal); %Turning received mapped signal into a sequence of bits
    numberOfErrors = 0;

    %Comparing original bits and recovered bits
    for n = 1:length(bits)
       if bits(n)~=recovBits(n)
           numberOfErrors = numberOfErrors + 1;
       end
    end
    
    ebN0_dec = 10^(ebN0_dB/10); %Decimal Eb/N0
    i = ebN0_dB + 1; %Indexing for the plot

    BER(i) = numberOfErrors / length(bits); %BER for current Eb/N0
    TBER(i) = qfunc(sqrt(2 * ebN0_dec)); %Theoretical BER

    %Plots
    set(empiricalPlot, 'XData', ebN0_range(1:length(BER)), 'YData', BER);
    set(theoreticalPlot, 'XData', ebN0_range(1:length(TBER)), 'YData', TBER);
    pause(0.1);
end
    %Plot setup
hold off;
legend('Empirical', 'Theoretical');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
grid on;
title('BER vs. E_b/N_0 for QPSK');


% This function will convert our Nx2 matrix of symbols into Nx1 vector of gray mapped symbols
function grayMapped = grayMap(bits)
    rows = size(bits, 1);
    grayMapped = zeros(rows, 1);
    for n = 1:rows
        bitPair = bits(n, :);
        if isequal(bitPair, [0, 0])
            grayMapped(n) = 1+1i;
        elseif isequal(bitPair, [1, 0])
            grayMapped(n) = -1+1i;
        elseif isequal(bitPair, [1, 1])
            grayMapped(n) = -1-1i;
        elseif isequal(bitPair, [0, 1])
            grayMapped(n) = 1-1i;
        end
    end
end

% This function will be called multiple times for each dB value of Eb/N0 and return signal injected with noise
function n = getNoiseAndSignal(ebN0_dB, noise, signal)
    ebN0_dec = 10^(ebN0_dB/10); %Converted to decimal
    noise_scaling = 1 / sqrt(2*ebN0_dec); %obtained scalar
    totalAWGN = noise_scaling.*noise; %scaled AWGN
    n = signal + totalAWGN; %injecting noise into our signal
end

% This function will recover the transmitted signal as it was mapped before AWGN
function reM = remap(sig)
    size = length(sig);
    reM = zeros(size, 1); 
    %since we mapped our signal as +-1+-1i, we can just compare the signs
    for n = 1:size
        r = real(sig(n));
        i = imag(sig(n));
        if (r>0 && i>0)
            reM(n) = 1+1i;
        elseif (r<0 && i>0)
            reM(n) = -1+1i;
        elseif (r<0 && i<0)
            reM(n) = -1-1i;
        elseif (r>0 && i<0)
            reM(n) = 1-1i;
        end
    end
end

% This function will demodulate the recovered signal and output bits
function deM = demodulate(sig)
    %Making Nx2 matrix that will store real component's bits in 1st column and the imaginary component's bit in the 2nd column
    N = length(sig);
    deM1 = zeros(length(sig), 2);
    for n = 1:N
        r = real(sig(n));
        i = imag(sig(n));

        if r < 0
            deM1(n, 1) = 1;
        else
            deM1(n, 1) = 0;
        end
        if i < 0
            deM1(n, 2) = 1;
        else
            deM1(n, 2) = 0;
        end
    end
    %Turning Nx2 matrix into a row vector
    deM = reshape(deM1', 1, []);
end

%Q function
function q = qfunc(x)
    q = 0.5 * erfc(x / sqrt(2));
end
