
% Function to define a magnet to work with calcMagField

% James O'Connell 17th December 2019

function mag = createMagnet(varargin)

for i = 1:nargin-1
    
    if strcmp(varargin{i},'vertices')
        mag.vertices = varargin{i+1};
    end
    
    if strcmp(varargin{i},'mur')
        mag.mur = varargin{i+1};
    end
    
    if strcmp(varargin{i},'magnetisation')
        mag.magnetisation = varargin{i+1};
    end



end