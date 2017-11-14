function Q = PrintContents(FN)

if exist('fullfilename', 'file'), FFN = fullfilename(FN, cd);
else, FFN = FN;
end

fNames = fieldnames(Q);

for n = 1:length(fNames)
    disp(Q.(fNames{n}))
end