function [ accuracy ] = test( Mapping, Testing_Data, Type )
% Function to test accuracy of the trained neural network according
% to the testing data
% Mapping is the CMAC Mapping created using the create function. It consists of the
% input vector, the look up table, the corrected weights, and the number of
% association cells linked to each input vector.
% Testing_Data is the data to be used to test the trained CMAC.

% tic

if isempty(Mapping) || isempty(Testing_Data)
    accuracy = NaN;
    return
end

% Input w.r.t. Input vectors
input  = zeros(length(Testing_Data),2);
for i=1:length(Testing_Data)
    if Testing_Data(i,1) > Mapping{1}(end)
        input(i,1) = length(Mapping{1});
    elseif Testing_Data(i,1) < Mapping{1}(1)
        input(i,1) = 1;
    else
        temp = (length(Mapping{1})-1)*(Testing_Data(i,1)-Mapping{1}(1))/(Mapping{1}(end)-Mapping{1}(1)) + 1;
        input(i,1) = floor(temp);
        if (ceil(temp) ~= floor(temp)) && Type
            input(i,2) = ceil(temp);
        end
    end
end

% Accuracy
numerator = 0;
denominator = 0;
for i=1:length(input)
    if input(i,2) == 0
        output = sum(Mapping{3}(find(Mapping{2}(input(i,1),:))));
        numerator = numerator + abs(Testing_Data(i,2)-output);
        denominator = denominator + Testing_Data(i,2) + output;
    else
        d1 = norm(Mapping{1}(input(i,1))-Testing_Data(i,1));
        d2 = norm(Mapping{2}(input(i,2))-Testing_Data(i,1));
        output = (d2/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,1),:))))...
               + (d1/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,2),:))));
        numerator = numerator + abs(Testing_Data(i,2)-output);
        denominator = denominator + Testing_Data(i,2) + output;
    end
    B(i) = output;
end
error = abs(numerator/denominator);
accuracy = 100 - error;

[A,I] = sort(Testing_Data(:,1));
B = B(I);
plot(A,B);

% toc

end