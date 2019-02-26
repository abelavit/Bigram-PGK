
% This code was developed using R2016a Matlab version (windows platform)
% Compare the samples in test dataset with already obtained labels from Phogly–PseAAC webserver at http://app.aporc.org/Phogly-PseAAC/
% The result of this code is stored in the variable Results_Avg_Phogly_PseAAC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

load Phogly_PseAAC_Result % Contains the prediction by the Phogly–PseAAC predictor
load test

for k = 1:10

    st = num2str(k);
    test_name = strcat('Fold', st, '_10');
    fold_file = [eval(test_name)];
    
    FN = 0;
    FP = 0;
    TN = 0;
    TP = 0;
    
    for i=1:size(fold_file,1)

        phosphoglycerylation_match = 0; % Put a value of zero and only put it to one when match has occured below in the code  
        
        for j=1:size(result_file,1)
       
           if (strcmp(fold_file{i,1}, result_file{j,2}) == true) && (result_file{j,3} == fold_file{i,4})

                phosphoglycerylation_match = 1; % Match has been found

           end
           
           yes_no{i,1} = phosphoglycerylation_match; % store the predicted label

        end
        
        % Obtaining the FN, FP, TN and TP values:
        
        if (fold_file{i,3} == '1') && (phosphoglycerylation_match == 1) % Actual and predicted are 1 

            TP = TP + 1;
            
        elseif (fold_file{i,3} == '1') && (phosphoglycerylation_match == 0) % Actual is 1 and predicted is 0

            FN = FN + 1;

        elseif (fold_file{i,3} == '0') && (phosphoglycerylation_match == 1) % Actual is 0 and predicted is 1

            FP = FP + 1;

        elseif (fold_file{i,3} == '0') && (phosphoglycerylation_match == 0) % Actual and predicted are 0

            TN = TN + 1;

        end

    end
    sen = TP/(TP+FN);
    spe = TN/(TN+FP);
    pre = TP/(TP+FP);
    accuracy = (TN+TP)/(FN+FP+TN+TP);
    mcc = ((TN*TP)-(FN*FP))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    Results_Phogly_PseAAC(1,k) = sen; 
    Results_Phogly_PseAAC(2,k) = spe; 
    Results_Phogly_PseAAC(3,k) = pre; 
    Results_Phogly_PseAAC(4,k) = accuracy; 
    Results_Phogly_PseAAC(5,k) = mcc;
    
    confusion_matrix(1,k) = FN;
    confusion_matrix(2,k) = FP;
    confusion_matrix(3,k) = TN;
    confusion_matrix(4,k) = TP;
    
    fold_label = cell2mat(fold_file(:,3)); % Save true labels (for AUC calculation)
    Test_label = str2num(fold_label); % converting type
    fold_true_label{k} = Test_label; % converting type
    
    fold_pred_label{k} = cell2mat(yes_no); % Save predicted labels (for AUC calculation)
    
    clear yes_no
    
end
metrics = {'sensitivity'; 'specificity'; 'precision'; 'accuracy'; 'mcc'};
Results_Avg_Phogly_PseAAC = [metrics, num2cell(sum(Results_Phogly_PseAAC,2)/10)];

% Display 10 fold CV result of Phogly_PseAAC
Ten_Fold_CV = [metrics, Results_Avg_Phogly_PseAAC(:,2)]

% AUC calculation 
classes = [fold_true_label{1}; fold_true_label{2}; fold_true_label{3}; fold_true_label{4}; fold_true_label{5}; fold_true_label{6}; fold_true_label{7}; fold_true_label{8}; fold_true_label{9}; fold_true_label{10}];
scores = [fold_pred_label{1}; fold_pred_label{2}; fold_pred_label{3}; fold_pred_label{4}; fold_pred_label{5}; fold_pred_label{6}; fold_pred_label{7}; fold_pred_label{8}; fold_pred_label{9}; fold_pred_label{10}];
[X,Y,T,AUC] = perfcurve(classes,scores,1);
AUC

%{
% Write to text file the actual label and predictions to calculate AUC

filename = fullfile('C:\Users\User\Desktop\AUC', 'AUC_Data_Phogly_PseAAC.txt'); % Change directory here
fid = fopen(filename,'w');
for l=1:size(AUC_Data_Phogly_PseAAC, 1)
    The_String = strcat(AUC_Data_Phogly_PseAAC{l,1}, ',', AUC_Data_Phogly_PseAAC{l,2}, ',', AUC_Data_Phogly_PseAAC{l,3});
    fprintf(fid, '%s\n', The_String);
end
fclose(fid)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
