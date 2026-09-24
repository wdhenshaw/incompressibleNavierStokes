function font( varargin )
% valid fontsize from 3 and up
% valid styles 0, 1, 2 for plain, bold, and italic. (string 'bold' ect.
% also ok.
% valid fonts 'Monospaced', 'SansSerif', 'Serif'
% wdh: Preferences->Fonts-> Desktop text font : use system font
%      font(13,'Lucida Grande');

nargin = length(varargin);
defaultText = 'Monospaced';
defaultSize = 14;
PLAIN = 0;
BOLD = 1;
ITALIC = 2;
defaultWeight = 0;
validWeights = {'plain','bold','italic'};
checkWeights = @(x) any(strcmp(x,validWeights));
validTexts = {'Monospaced', 'SansSerif','Serif'};
checkTexts = @(x) any(strcmp(x,validTexts));
if nargin == 0
  weight  = defaultWeight;
  text = defaultText;
  size = defaultSize;
  com.mathworks.services.FontPrefs.setCodeFont(...
    java.awt.Font(text,weight,size));
  return
end
set = [];
for i = 1:nargin
  arg = varargin{i};
  if (isnumeric(arg) || ~isempty(str2num(arg)) )
    if ~isnumeric(arg)
      arg = str2num(arg);
    end
    if arg <= 2 && (any(set~=2) || isempty(set))
      weight = arg;
      set = [set,2];
    else
      size = arg;
      set = [set,3];
    end
  else
    if checkTexts(arg)
      text = arg;
      set = [set, 1];
    elseif checkWeights(arg)
      for i = 1:3
        if strcmp(arg,validWeights{i})
          weight = i-1;
          set = [set,2];
        end
      end
    end
  end
end

if ~any(1==set)
  text = defaultText;
end
if ~any(2==set)
  weight = defaultWeight;
end
if ~any(3==set)
  size = defaultSize;
end
com.mathworks.services.FontPrefs.setCodeFont(...
    java.awt.Font(text,weight,size));
return
end
