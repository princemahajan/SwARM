function [Clm,Slm] = DenormCS(CSfilename, l_max)
%DenoemCS Denormalize C and S coefficients and return them as arrays
%
% Format of File: '%6s   %f %f   %21f%21f'
%
% Author: Bharat Mahajan (https://github.com/princemahajan)


% read and ignore the lines preceded by #
fid = fopen(CSfilename);
if fid == -1
    disp(strcat('DenormCS: check coeff file path ', CSfilename));
end

ctr = 0;
while true
    % read line
    line = fgets(fid);
    
    if line(1) ~= '#'
        % move the counter back to the begining of the line
        fseek(fid, -length(line), 'cof');
        break;
    end
    ctr = ctr + 1;
    
end

% read from coeff file
CSlm = textscan(fid,'%6s   %f %f   %f %f');
fclose(fid);

% memory allocation
Clm = zeros(l_max+1,l_max+1);
Slm = Clm;

% Normalization constant
pilmk = @(l,m,k) sqrt(factorial(l + m)/(factorial(l - m)*k*(2*l + 1)));

% start from degree 2
lctr0 = 2;

% Denormalize
for ctr = 1:length(CSlm{2})


    % degree ctr
    l = CSlm{2}(ctr);
    lctr = lctr0 + (l - 1);

    if l > l_max || isnan(l) == true
        break;
    end
    
    % order counter
    m = CSlm{3}(ctr);
    mctr = CSlm{3}(ctr) + 1;
    
    k = 2;
    if m == 0
        k = 1;
    end
    
    Clm(lctr,mctr) = CSlm{4}(ctr)/pilmk(l,m,k);
    Slm(lctr,mctr) = CSlm{5}(ctr)/pilmk(l,m,k);
    
%     disp(l);
%     disp(m);
    
end
    
% Make Slm zeros for m=0
Slm(:,1) = 0;

end

