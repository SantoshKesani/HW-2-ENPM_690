function [ Mapping, iteration, finalError, t ] = train( parentMapping, Training_Data, E, Type )
% Function to train CMAC according to the training data.
% Implies, Function gives the CMAC architecture with the corrected weights.
% Mapping is the CMAC Mapping created using the create function.
% Training_Data is the data to be used to train CMAC.
% E is the acceptable error. This error is in terms of the training data.
% Type = determines discrete or continuous

tic;

Mapping = parentMapping;
if isempty(Mapping) || isempty(Training_Data) || isempty(E)
    return
end

% Input w.r.t. Input vectors
input  = zeros(length(Training_Data),2);
for i=1:length(Training_Data)
    if Training_Data(i,1) > Mapping{1}(end)
        input(i,1) = length(Mapping{1});
    elseif Training_Data(i,1) < Mapping{1}(1)
        input(i,1) = 1;
    else
        temp = (length(Mapping{1})-1)*(Training_Data(i,1)-Mapping{1}(1))/(Mapping{1}(end)-Mapping{1}(1)) + 1;
        input(i,1) = floor(temp);
        if (ceil(temp) ~= floor(temp)) && Type
            input(i,2) = ceil(temp);
        end
    end
end

% compute output for each input and accordinglB adjust weights until the
% specified number of iterations is achieved
eta = 0.025; % learning rate
error = Inf;
iteration = 0;
count = 0;
while (error > E)&&(2*count <= iteration)
    old_err = error;
    iteration = iteration + 1;
    
    % compute output for each input and accordinglB adjust weights
    for i=1:length(input)
        if input(i,2) == 0
            output = sum(Mapping{3}(find(Mapping{2}(input(i,1),:))));
            error = eta*(Training_Data(i,2)-output)/Mapping{4};
            Mapping{3}(find(Mapping{2}(input(i,1),:))) = Mapping{3}(find(Mapping{2}(input(i,1),:))) + error;
        else
            d1 = norm(Mapping{1}(input(i,1))-Training_Data(i,1));
            d2 = norm(Mapping{1}(input(i,2))-Training_Data(i,1));
            output = (d2/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,1),:))))...
                    + (d1/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,2),:))));
            error = eta*(Training_Data(i,2)-output)/Mapping{4};
            Mapping{3}(find(Mapping{2}(input(i,1),:))) = Mapping{3}(find(Mapping{2}(input(i,1),:)))...
                                                    + (d2/(d1+d2))*error;
            Mapping{3}(find(Mapping{2}(input(i,2),:))) = Mapping{3}(find(Mapping{2}(input(i,2),:)))...
                                                    + (d1/(d1+d2))*error;            
        end
    end

    % compute final error
    numerator = 0;
    denominator = 0;
    for i=1:length(input)
        if input(i,2) == 0
            output = sum(Mapping{3}(find(Mapping{2}(input(i,1),:))));
            numerator = numerator + abs(Training_Data(i,2)-output);
            denominator = denominator + Training_Data(i,2) + output;
        else
            d1 = norm(Mapping{1}(input(i,1))-Training_Data(i,1));
            d2 = norm(Mapping{1}(input(i,2))-Training_Data(i,1));
            output = (d2/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,1),:))))...
                   + (d1/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,2),:))));
            numerator = numerator + abs(Training_Data(i,2)-output);
            denominator = denominator + Training_Data(i,2) + output;
        end
    end
    error = abs(numerator/denominator);
    if abs(old_err - error) < 0.00001
        count = count + 1;
    else
        count = 0;
    end
end
iteration = iteration - count;

% compute final error
numerator = 0;
denominator = 0;
for i=1:length(input)
    if input(i,2) == 0
        B(i) = sum(Mapping{3}(find(Mapping{2}(input(i,1),:))));
        numerator = numerator + abs(Training_Data(i,2)-B(i));
        denominator = denominator + Training_Data(i,2) + B(i);
    else
        d1 = norm(Mapping{1}(input(i,1))-Training_Data(i,1));
        d2 = norm(Mapping{1}(input(i,2))-Training_Data(i,1));
        B(i) = (d2/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,1),:))))...
               + (d1/(d1+d2))*sum(Mapping{3}(find(Mapping{2}(input(i,2),:))));
        numerator = numerator + abs(Training_Data(i,2)-B(i));
        denominator = denominator + Training_Data(i,2) + B(i);
    end
end
finalError = abs(numerator/denominator);
[A,I] = sort(Training_Data(:,1));
B = B(I);
% plot(A,B);

t = toc;

end