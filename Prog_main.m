% clear all
close all
clc

A = (linspace(1,25))';
B = sin(A);
% Training_Data_3 = [A(1:70),B2(1:70)];
% Testing_Data_3 = [A(71:100),B2(71:100)];
% CMAC_Mapping = create_2(A,35);
% [CMAC_Mapping_3,iteration_3,error_3,t_3] = train_2(CMAC_Mapping,Training_Data_3,1);
% accuracy_3 = test_2(CMAC_Mapping_3,Testing_Data_3);
% iteration = zeros(2,34);
% iter = zeros(2,34);
% accuracy = zeros(2,34);
% acc = zeros(2,34);
% t = zeros(2,34);
% T = zeros(2,34);
% for j=1:100
    I = randperm(100);
    Training_Data = [A(I(1:70)),B(I(1:70))];
    Testing_Data = [A(I(71:100)),B(I(71:100))];
    for i=1:34
        CMAC_Mapping = create(A,35,i);
        figure
        plot(A,B);
        hold on
        [Mapping,iter(1,i),~,T(1,i)] = train(CMAC_Mapping,Training_Data,0,0);
        acc(1,i) = test(Mapping,Testing_Data,0);
        hold off
        legend('Main Function','Test Output');
        title(['Overlap = ' num2str(i)]);
        figure
        plot(A,B);
        hold on
        [Mapping,iter(2,i),~,T(2,i)] = train(CMAC_Mapping,Training_Data,0,1);
        acc(2,i) = test(Mapping,Testing_Data,1);
        hold off
        legend('Original Curve','Testing DataB');
        title(['CellNo = ' num2str(i)]);
    end
% iteration = iteration + iter;
% accuracy = accuracy + acc;
% t = t + T;
% end
% iteration = iteration/j;
% accuracy = accuracy/j;
% t = t/j;