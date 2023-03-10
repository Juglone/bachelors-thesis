% The code is based on the examples found in
% "Mmwave Radar Device ADC Raw Data Capture", download at:
% https://www.ti.com/lit/an/swra581b/swra581b.pdf

% retVal will contain a matrix, row 1 contains all of the data from the first receiver, row 2 from the second receiver, row 3
% from the third receiver, and row 4 from the fourth receiver. In cases where certain receivers are disabled,
% the corresponding rows will be removed. Each row will contain a number of columns equal to the number
% of ADC samples per chirp multiplied by the total number of chirps. The columns are organized by chirps.
% For example, if there are 256 ADC samples and a total of 8 chirps, each row will contain 2048 columns.
% The first 256 columns will correspond to the first chirp, the next 256 columns to the second chirp, and so on. 
function [retVal] = readDCA1000(fileName,samples_per_chirp)

    %% global variables
    % change based on sensor config
    numADCSamples = samples_per_chirp; % number of ADC samples per chirp
    numADCBits = 16; % number of ADC bits per sample
    numRX = 1; % number of receivers
    numLanes = 2; % do not change. number of lanes is always 2
    isReal = 0; % set to 1 if real only data, 0 if complex data0

    %% read file
    % read .bin file
    fid = fopen(fileName,'r');
    adcData = fread(fid, 'int16');

    % if 12 or 14 bits ADC per sample compensate for sign extension
    if numADCBits ~= 16
        l_max = 2^(numADCBits-1)-1;
        adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
    end
    fclose(fid);
    fileSize = size(adcData, 1);

    % real data reshape, filesize = numADCSamples*numChirps
    if isReal
        numChirps = fileSize/numADCSamples/numRX;
        LVDS = zeros(1, fileSize);
        %create column for each chirp
        LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
        %each row is data from one chirp
        LVDS = LVDS.';
    else
        % for complex data
        % filesize = 2 * numADCSamples*numChirps
        numChirps = fileSize/2/numADCSamples/numRX;
        LVDS = zeros(1, fileSize/2);
        %combine real and imaginary part into complex data
        %read in file: 2I is followed by 2Q
        counter = 1;
        for i=1:4:fileSize-1
            LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
            LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
            counter = counter + 2;
        end
        % create column for each chirp
        LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
        %each row is data from one chirp
        LVDS = LVDS.';
    end

    %organize data per RX
    adcData = zeros(numRX,numChirps*numADCSamples);
    for row = 1:numRX % kan egentligen strunta i den här då vi bara har en mottagarantenn
        for i = 1: numChirps
            adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = ...
                LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);   
        end
    end

    % return receiver data
    retVal = adcData;
end
