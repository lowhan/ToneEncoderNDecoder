%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      Assignment: Tone encoder and decoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frequencies = containers.Map;                           % Assigned frequencies for each digit
frequencies('0') = [320, 750];
frequencies('1') = [410, 875];
frequencies('2') = [515, 912];
frequencies('3') = [640, 1080];
fs = 8000;                                              % sampling frequency
input = '1032';                                         % input stream of digits
toneduration = 0.5;                                     % 0.5 seconds for each tone                            
filter_threshold = 1200;                                % cutoff frequency on decode

% conversion of continuous time domain to discrete time domain
t = linspace(0, toneduration, toneduration*fs);         % continuous time domain
% (8000 occurrences per 1 sec) * 0.5s = 4000 occurrences in 0.5 duration = 0.000125 second each points
% t = nTs > n = t/Ts > n = t*fs
% n = (1:(toneduration*fs));  
n = t*fs;                                               % discrete time domain for 0.5s
number_figure = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tone encoder: Dual Tone Generator implementation
encoded_tone = [];
for i = 1:length(input)                     % Iterate over each digit in the input stream '1032'
    digit = input(i);
    freq = frequencies(digit);              % load the digits's frequencies pair
    f1 = freq(1);
    f2 = freq(2);
    w1 = 2*pi*f1/fs;
    w2 = 2*pi*f2/fs;
    % Generate the mixed tone for the digit (0 - 0.5s @ 1 - 4000) as x(n)
    tone1 = sin(w1*n);              
    tone2 = sin(w2*n);
    output = tone1 + tone2;                 % Combine tone 1 and tone 2
    encoded_tone = [encoded_tone, output];  % Produce output

    % Plot the graph of encoded tone pair
    number_figure = number_figure + 1;
    figure('Name','Encoding the tone pairs');
    tiledlayout(3,1);
    nexttile(1);
    plot(n,tone1);
    xlabel('n'); ylabel('Amplitude, A');
    title("Tone 1, f1");
    nexttile(2);
    plot(n,tone2);
    xlabel('n'); ylabel('Amplitude, A');
    title("Tone 2, f2");
    nexttile(3);
    plot(n,output);
    xlabel('n'); ylabel('Amplitude, A');
    title("Combined Output");
    sgtitle("Encoding of digit " + digit);
    figure(number_figure);
end

% Play the encoded tone pairs
sound(encoded_tone, fs);

% 4000 occurrences * 4 digits = 16000 occurrences
% 4 digits each take 0.5s on tone duration = 2s used
% same to (16000 occurrences / (8000 occurrences per 1 sec)) = 2s

% Range of discrete time for 2s
N = (1:length(encoded_tone)); 

% Plot the graph of the mixed signal of input stream '1032'
number_figure = number_figure + 1;
figure('Name','Tone Encoder');
plot(N,encoded_tone);
xlabel('n'); ylabel('Amplitude, A');
title("Encoded Dual Tone Multi Frequency Signal, x(n)");
figure(number_figure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Tone decoder using Goertzel algorithm
decoded_digits = [];
num_pairs = length(encoded_tone)/(toneduration*fs);
Goertzel_range = [320 410 515 640 750 875 912 1080];        % limit the range for performing Goertzel algorithm
Goertzel_index = round(Goertzel_range/fs*max(n)) + 1;       % Index for each frequency
Goertzel_threshold = 0.1;

number_figure = number_figure + 1;
figure('Name','Tone Decoder using Goertzel algorithm');
figure(number_figure);
for i = 1:num_pairs
    % Separate first 4000 points of a dual tone pair (1 ~ 4000, 4001 ~ 8000, 8001 ~ 12000, 12001 ~ 16000)
    tone_segment = encoded_tone((i-1)*(toneduration*fs)+1 : i*(toneduration*fs));
    
    % Compute DFT of a tone segment using goertzel algorithm
    tic;
    spectrum = abs(goertzel(tone_segment/fs,Goertzel_index));    % amplitude is normalized and absolute

    % Find the first two dominant frequencies
    index = find(spectrum >= Goertzel_threshold);                % find the first two peak that has atleast 0.1 amplitude and above
    first_peak = Goertzel_range(index(1));        
    second_peak = Goertzel_range(index(2));   

    % Comparing the correspond first_peak and second_peak with the samples
    sample = keys(frequencies);
    for digit = 1:length(sample)
        target = frequencies(sample{digit});
        difference = min(norm(target(1) - first_peak), norm(target(2)-second_peak));
        if difference == 0
            closest_digit = sample{digit};
            min_difference = difference;
        end
    end
    % Add the decoded digit into stream
    decoded_digits = [decoded_digits, closest_digit]; 
    toc

    % Show the result of each digit in frequency domain
    nexttile(i);
    stem(Goertzel_range,spectrum);
    xlabel('frequency, k'); ylabel('Amplitude, A');
    title("Digit " + closest_digit);
    sgtitle("Decoded Dual Tone Multi Frequency Signal, x(k) using Goertzel Algorithm");
end

% Display the decoded stream in command window for using goertzel algorithm
disp("Decoded stream using goertzel algorithm: " + decoded_digits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           Tone decoder using DFT/FFT
decoded_digits_fft = [];   
num_pairs = length(encoded_tone)/(toneduration*fs);

number_figure = number_figure + 1;
figure('Name','Tone Decoder using DFT/FFT');
figure(number_figure);
for i = 1:num_pairs
    % Separate first 4000 points of a dual tone pair (1 ~ 4000, 4001 ~ 8000, 8001 ~ 12000, 12001 ~ 16000)
    tone_segment = encoded_tone((i-1)*(toneduration*fs)+1 : i*(toneduration*fs));
    
    % Compute DFT of a tone segment using fft
    tic;
    scope = abs(fft(tone_segment/fs));         % amplitude is normalized and absolute
    scope_length = length(scope);              % length
    scope_length_range = (1:scope_length);     % range
    scope_resolution = fs/scope_length;        % ratio between scope and fs
    [~, B] = sort(scope,'ascend');             % rearrange the output of separated dual tone pairs in ascending order
    filterB = B.*(B<filter_threshold);         % cancel the frequencies that above the threshold
    indexB = (nonzeros(filterB)).';            % remove the zeros            

    % Find the first two dominant frequencies
    first_peak = indexB(end);        
    second_peak = indexB(end-1);   

    % Since the resolution is half, multiply with the ratio to match with initial tone frequency
    first_frequency_tone = first_peak*scope_resolution;
    second_frequency_tone = second_peak*scope_resolution;
    
    % Comparing the correspond first_peak and second_peak with the samples
    min_difference = inf;
    sample = keys(frequencies);
    for digit = 1:length(sample)
        target = frequencies(sample{digit});
        difference = min(norm(target(1) - second_frequency_tone), norm(target(1)- first_frequency_tone));
        if difference < min_difference 
            closest_digit = sample{digit};
            min_difference = difference;
        end
    end
    % Add the decoded digit into stream
    decoded_digits_fft = [decoded_digits_fft, closest_digit]; 
    toc

    % Show the result of each digit in frequency domain
    nexttile(i);
    plot(scope_length_range,scope);
    xlabel('frequency, k'); ylabel('Amplitude, A');
    title("Digit " + closest_digit);
    sgtitle("Decoded Dual Tone Multi Frequency Signal, x(k) using DFT/FFT");
end

% Display the decoded stream in command window for using DFT/FFT
disp("Decoded stream using DFT/FFT: " + decoded_digits_fft);
