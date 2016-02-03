function wienerModel(fnum,varargin)

path(path,'../analyze');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manage options
SHOW_PLOT = 0;
SHOW_R2 = 0;
SAVE_RESULTS = 0;
for i = 1:nargin-1
    switch varargin{i}
        case 'plot'
            SHOW_PLOT = 1;
        case 'r2'
            SHOW_R2 = 1;
        case 'save'
            SAVE_RESULTS = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load input data
fNameIn = sprintf('../results/d%05iinp.bin',fnum);
fid = fopen(fNameIn,'rb');
fread(fid,50,'char'); % dataType
fread(fid,1,'double');% tMax
fread(fid,5,'double'); %dataParam
nNeurons = fread(fid,1,'double');
dt = fread(fid,1,'double');
nSamples = fread(fid,1,'double');
posn = fread(fid,nSamples,'double');
N = fread(fid,[nNeurons,nSamples],'char');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the analysis
pctTrain = 0.5;
nLags = 20;
windowDur = 50e-3; % seconds
winWidth  = windowDur / dt; % samples

nFRblocks = floor(nSamples/winWidth);
FR = zeros(nNeurons,nFRblocks);
for i = 1:nNeurons
    FR(i,:) = sum( reshape(N(i,:),winWidth,[]) , 1 );
end
FR = FR';

rTrain = 1 : floor(pctTrain * nFRblocks);
rTest  = rTrain(end)+1 : nFRblocks;

frTrain = buildLagMatrix(FR(rTrain,:),nLags);
frTest  = buildLagMatrix(FR(rTest ,:),nLags);

posn2 = mean(reshape(posn,winWidth,[]),1);

posnTrain = reshape(posn2(rTrain(nLags:end)),[],1);
posnTest  = reshape(posn2(rTest (nLags:end)),[],1);

X = frTrain; Y = posnTrain;
% A = X\Y;
% % A = inv(X' X)X' Y
A = pinv(X' * X) * X' * Y;

X = frTest;
Y = X*A;
posnEstd = Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SHOW_PLOT
    clf;
    t = (1:length(posnEstd)) * dt * winWidth;
    plot(t,posnEstd,'tag','not_thick'); hold on;
    plot(t,posnTest,'r');
    xtight; xlabel('time(s)'); legend('estd','true');
    bigText;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2 = computeR2(posnTest,posnEstd);
if SHOW_R2
    disp(r2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SAVE_RESULTS
    fNameOut = sprintf('../results/d%05iwnr.bin',fnum);
    fid = fopen(fNameOut,'wb');
    fwrite(fid,sprintf('%-10s','linVer1'),'char');
    fwrite(fid,sprintf('%-50s',datestr(clock)),'char');
    fwrite(fid,sprintf('%-50s',fNameIn),'char');
    fwrite(fid,pctTrain,'double');
    fwrite(fid,nLags,'double');
    fwrite(fid,windowDur,'double');
    fwrite(fid,length(posnEstd),'double');
    fwrite(fid,dt,'double');
    fwrite(fid,r2,'double');
    fwrite(fid,posnEstd,'double');
    fwrite(fid,posnTest,'double');
    fclose(fid);
end
