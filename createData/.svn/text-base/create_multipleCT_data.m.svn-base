function inputData = create_multipleCT_data(nNeurons,T,varargin)

FLAT_MU = 0;
if nargin>2
    for i = 3:nargin
        ind = i -2;
        if strcmp(varargin{ind},'flat')
            FLAT_MU = 1;
        end
    end
end

dt = 2e-3; %seconds
t = dt:dt:T; %seconds

th_p = repmat( 2*pi*rand(nNeurons,1)  , 1 , length(t) );
fr_L = repmat( 30*rand(nNeurons,1)    , 1 , length(t) );
fr_H = repmat( 30*rand(nNeurons,1)+70 , 1 , length(t) );

model = @(th_p,fr_L,fr_H,x) ((fr_H - fr_L)/2).*cos(x - th_p) + (fr_H + fr_L)/2;

posn = (pi/2)*sin(2*pi/5 * t) + pi;
lambda = model(th_p,fr_L,fr_H,repmat(posn,nNeurons,1));
N = poissrnd(lambda*dt);
N = N>0;

inputData.dt = dt;
inputData.t = t;
inputData.posn = posn;
inputData.N = N;
inputData.x.th_p = th_p;
inputData.x.fr_L = fr_L;
inputData.x.fr_H = fr_H;
inputData.nNeurons = nNeurons;

% save data02 inputData;
