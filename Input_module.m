function [ Mapping ] = create( A, noWts, CellNo )
% Function to link input to the association cells & assigning weights to 1
% A = input values
% noWts = No. of weights
% CellNo = No. of Association cells
% Note that, the CellNo option implies that overlap between the
% consecutive input vectors is equal to CellNo-1.

if (CellNo > noWts) || (CellNo < 1) || (isempty(A))
    Mapping = [];
    return
end

% Input Vector
A = linspace(min(A),max(A),noWts-CellNo+1)';

% Look up table
LUtable = zeros(length(A),noWts);
for i=1:length(A)
    LUtable(i,i:CellNo+i-1) = 1;
end

% Wts.
W = ones(noWts,1);

% Mappingping
Mapping = cell(3,1);
Mapping{1} = A;
Mapping{2} = LUtable;
Mapping{3} = W;
Mapping{4} = CellNo;

end