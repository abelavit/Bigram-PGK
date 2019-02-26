
% This code was developed using R2016a Matlab version (windows platform)
% Compare the samples in test dataset with already obtained labels from CKSAAP_PhoglySite matlab software available at https://github.com/juzhe1120/Matlab_Software/blob/master/CKSAAP_PhoglySite_Matlab_Software.zip
% The result of this code is stored in the variable Results_Avg_CKSAAP_PhoglySite 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

load CKSAAP_PhoglySite_Result % Contains the prediction by the CKSAAP_PhoglySite predictor
load test

Result_File = aaa(3:size(aaa,1),:); % We are not intereted in titles of columns hence rows starting at 3

for k = 1:10

    st = num2str(k);
    test_name = strcat('Fold', st, '_10');
    fold_file = [eval(test_name)];
    
    FN = 0;
    FP = 0;
    TN = 0;
    TP = 0;

    for i=1:size(fold_file,1)

        for j=1:size(Result_File,1)

            if (strcmp(Result_File(j,1), fold_file{i,1}) == true ) && (Result_File{j,2} == fold_file{i,4})

                if Result_File{j,3} == -1
                    
                    yes_no{i,1} = 0; % store label of 0
                    
                else
                    
                    yes_no{i,1} = Result_File{j,3}; % store label of 1
                    
                end
            end
        end

        % Obtaining the FN, FP, TN and TP values:
        
        if yes_no{i,1} == 1 % Predicted is 1

            if fold_file{i,3} == '1' % Actual is 1

                TP = TP + 1;

            elseif fold_file{i,3} == '0' % Actual is 0

                FP = FP + 1;

            end

        elseif yes_no{i,1} == 0 % Predicted is 0

            if fold_file{i,3} == '1' % Actual is 1

                 FN = FN + 1;

            elseif fold_file{i,3} == '0' % Actual is 0

                 TN = TN + 1;

            end 
        end
    end
    sen2 = TP/(TP+FN);
    spe2 = TN/(TN+FP);
    pre2 = TP/(TP+FP);
    accuracy2 = (TN+TP)/(FN+FP+TN+TP);
    mcc2 = ((TN*TP)-(FN*FP))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    Results_CKSAAP_PhoglySite(1,k) = sen2; 
    Results_CKSAAP_PhoglySite(2,k) = spe2; 
    Results_CKSAAP_PhoglySite(3,k) = pre2; 
    Results_CKSAAP_PhoglySite(4,k) = accuracy2; 
    Results_CKSAAP_PhoglySite(5,k) = mcc2;
    
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
metrics = {'sensitivity'; 'specificity'; 'Precision'; 'accuracy'; 'mcc'};
Results_Avg_CKSAAP_PhoglySite = [metrics, num2cell(sum(Results_CKSAAP_PhoglySite,2)/10)];

% Display 10 fold CV result of CKSAAP_PhoglySite
Ten_Fold_CV = [metrics, Results_Avg_CKSAAP_PhoglySite(:,2)]

% AUC Calculation:

classes = [fold_true_label{1}; fold_true_label{2}; fold_true_label{3}; fold_true_label{4}; fold_true_label{5}; fold_true_label{6}; fold_true_label{7}; fold_true_label{8}; fold_true_label{9}; fold_true_label{10}];
scores = [fold_pred_label{1}; fold_pred_label{2}; fold_pred_label{3}; fold_pred_label{4}; fold_pred_label{5}; fold_pred_label{6}; fold_pred_label{7}; fold_pred_label{8}; fold_pred_label{9}; fold_pred_label{10}];
[X,Y,T,AUC] = perfcurve(classes,scores,1);
AUC
