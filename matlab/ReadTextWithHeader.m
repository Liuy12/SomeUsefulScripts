function data = ReadTextWithHeader( fileName, flag )
% read data into columns with header information.
%   data = ReadTextWithHeader( fileName, flag )
%   data = ReadTextWithHeader( fileName )
%
%   Output.
%       data.c:     column contents (n columns in cells)
%       data.colHeaders:    column header text.
%
%     if flag == 0, with header
%     if flag == 1, convert header to structure.
%       each column will be a field in the "data" structure. So make sure
%       the name of the column headers is reasonable...
%     if flag === 2
%       each column will be read without header.
%     if flag == 3
%       it will take first column as ID column, first row as name, and rest
%       will be converted to number.

% Yidong Chen.

if( nargin < 2 ), flag = 0; end;

fid = fopen( fileName, 'r' );
bLine = fgetl( fid );

expr = '([^\t]+).*?';
tok = regexp( bLine, expr, 'tokens' );
for i = 1:length(tok), colHeaders(i) = tok{i}; end;

% fix colHeader.
for i = 1:length( colHeaders )
    h = colHeaders{i};
    
    % space replaced with '_';
    if(flag~=3)
        idx = find( h == ' ' | h == '(' | h == '/' | h == '-' | h == '+' );
        if( ~isempty(idx) ), h(idx) = '_'; end
    end
    % '?' removed;
    idx = find( h == '?' | h == ')' );
    if( ~isempty(idx) ), h(idx) = []; end
    
    % '#' replaced with 'N';
    idx = find( h == '#' );
    if( ~isempty(idx) ), h(idx) = 'N'; end
    
    colHeaders{i} = h;
end

if( flag == 2 )
    % without header, what we got is only good for the number of columns!!
    fclose( fid );
    fid = fopen( fileName, 'r' );
end

switch( flag )
    case 0,         % default.
        strFormat = repmat( '%s', 1, length(colHeaders));
        c = textscan( fid, strFormat, 'delimiter', '\t' );
        
        data.c = c;
        data.colHeaders = colHeaders;
    case 1,
        strFormat = repmat( '%s', 1, length(colHeaders));
        c = textscan( fid, strFormat, 'delimiter', '\t' );

        % get all annotation related.
        for i = 1:length(colHeaders)
            data.(colHeaders{i}) = c{i};
        end
    case 2,
        strFormat = repmat( '%s', 1, length(colHeaders));
        c = textscan( fid, strFormat, 'delimiter', '\t' );

        data.c = c;
    case 3,
        strFormat = ['%s' repmat( '%f', 1, length(colHeaders)-1)];
        c = textscan( fid, strFormat, 'delimiter', '\t' );

        data.(colHeaders{1}) = c{1};
        data.names = colHeaders(2:end);
        data.value = zeros( length(c{1}), length( data.names ) );
        for i = 1:length( data.names )
            data.value(:,i) = c{i+1};
        end
end
fclose( fid );

return;

