function trimmedData = trimmedIDft(dftData, trimmedDataLength)
    untrimmedData = length(dftData)*ifft(dftData);
    trimmedData = untrimmedData(1:trimmedDataLength);
end