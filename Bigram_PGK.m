
clear all
close all
 
load train % Load Bigram-PGK train data
load test % Load Bigram-PGK test data

cd 'Give libsvm directory here'
% Change directory to where libsvm is saved. e.g. cd 'C:\Users\User\Desktop\libsvm-weights-3.22\matlab'

for i = 1:10

    clear Test_Data Train_Data test_label train_label Test_label Train_label Test Train
    st = num2str(i);
    test_name = strcat('Fold', st, '_10');
    train_name = strcat('Train', st);
    
    % Prepare the train, test, train label and test label sets
    Test_Data = [eval(test_name)];
    Train_Data = [eval(train_name)];
    
    Test = cell2mat(Test_Data(:,2)); % Get the features
    Train = cell2mat(Train_Data(:,2)); % Get the features

    test_label = cell2mat(Test_Data(:,3)); % Get the label (originally in string)
    train_label = cell2mat(Train_Data(:,3)); % Get the label (originally in string)
    Test_label = str2num(test_label); % Change string label to num
    Train_label = str2num(train_label); % Change string label to num

    % Train LibSVM-weights. Empty square brackets denote we are not supplying weights
    model=svmtrain([],Train_label,Train,['-s 0 -t 1 -c 1 -g 1']);
    
    % Test the classifier
    [pred,acc,prob_values]=svmpredict(Test_label,Test,model);
    
    % Save the predicted labels and the true labels
    prediction{i} = pred;
    True_label{i} = Test_label;

    % Obtaining the FN, FP, TN and TP values
    FN = 0;
    FP = 0;
    TN = 0;
    TP = 0;
    for j = 1:size(Test_label,1)
        if Test_label(j) == 1
            if pred(j) == 1
                TP = TP + 1;
            else
                FN = FN + 1; 
            end
        else
            if pred(j) == 1
                FP = FP + 1;
            else
                TN = TN + 1; 
            end
        end
    end
    
    % Calculating the performance metrics
    sen = TP/(TP+FN);
    spe = TN/(TN+FP);
    pre = TP/(TP+FP);
    accuracy = (TN+TP)/(FN+FP+TN+TP);
    mcc = ((TN*TP)-(FN*FP))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    Results(1,i) = sen; 
    Results(2,i) = spe; 
    Results(3,i) = pre; 
    Results(4,i) = accuracy; 
    Results(5,i) = mcc;
end
Results_Avg = sum(Results,2)/10; % Average the result to get 10 fold CV  

% AUC calculation: 

classes = [True_label{1}; True_label{2}; True_label{3}; True_label{4}; True_label{5}; True_label{6}; True_label{7}; True_label{8}; True_label{9}; True_label{10}];
scores = [prediction{1}; prediction{2}; prediction{3}; prediction{4}; prediction{5}; prediction{6}; prediction{7}; prediction{8}; prediction{9}; prediction{10}];

[X,Y,T,AUC] = perfcurve(classes,scores,1);

% Display the result
metrics = {'sensitivity'; 'specificity'; 'precision'; 'accuracy'; 'mcc'};
Ten_Fold_CV = [metrics, num2cell(Results_Avg)] 
AUC

