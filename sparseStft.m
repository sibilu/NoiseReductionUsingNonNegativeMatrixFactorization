function [sStft, freqVector, timeVector] = ...
        sparseStft(data, window, nOverlap, nDft, samplingFreq, lambda)
    nData = length(data);
    segmentLength = length(window);
    segmentTime = segmentLength/samplingFreq;
    nShift = segmentLength-nOverlap;
    nShiftTime = nShift/samplingFreq;
    nSegments = floor((nData-nOverlap)/nShift);
    timeVector = segmentTime/2+(1:nSegments)'*nShiftTime-nShiftTime/2;
    freqVector = samplingFreq*(0:nDft-1)'/nDft;
    
    sStft = nan(length(freqVector)/2+1,length(timeVector));
    dftData = @(segment) fft(segment, nDft);
    iDftData = @(dftSegment) trimmedIDft(dftSegment, segmentLength);
    idx = 1:segmentLength;
    for ii = 1:nSegments
        if mod(ii,10) == 0
            disp(['Processing segment ', num2str(ii), ' of ', ...
                num2str(nSegments)]);
        end
        windowedData = window.*data(idx);
        bpDft = bpd_salsa(windowedData, iDftData, dftData, nDft, lambda, ...
            500, 50);
        sStft(:,ii) = bpDft(1:nDft/2+1);
        idx = idx+nShift;
    end
end
