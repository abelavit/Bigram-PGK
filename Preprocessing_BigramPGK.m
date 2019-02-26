
% Algorithm below was used to extract features and create train and test samples for Bigram-PGK predictor
% The train sets after code execution are labelled Train1, Train2, ...,Train10
% The test sets after code execution are labelled Fold1_10, Fold2_10, ..., Fold10_10
% The train and test sets generated each time will be different but performance would be similar

clear all
close all

%-------------------------------------------------------------
% Data Extraction 
%-------------------------------------------------------------

load Phosphoglycerylationstruct

Unprocessed_data = DB_Phosphoglycerylation;

Field = size(Unprocessed_data,2); % Columns of unprocessed data. It is num of protein sequences 

z=0;
for a = 1:Field
    z = z + size(strfind(Unprocessed_data(a).seq{1},'K'), 2); % Finding total number of -ve and +ve samples
end

Final_Data = cell(z,4);
k=0;  

for l = 1:Field
    K_locations_field{l} = strfind(Unprocessed_data(l).seq{1},'K'); % Saves the locations of K found in protein sequences
end

window = 32; % 32 upstream and 32 downstream
knn_neighbor = 112; % knn is actually 111 as it uses distance with itself in the process

for m = 1:Field
    for n=1:size(K_locations_field{m},2) % loop in all the K locations
        
        k = k+1; % Increment to next data saving location
        Final_Data{k,1} = Unprocessed_data(m).name; % Save protein name at 1st position
        
        if K_locations_field{m}(n) <= window % Location K close to N terminus
            matrix_a = Unprocessed_data(m).pssm_prob(1:K_locations_field{m}(n)+ window,:); % Matrix containing pssm_prob from start of protein seq till 32 upstream of K
            matrix_b = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+ window,:); % Matrix containing pssm_prob from location 2x the location of K till end of upstream 
            matrix_c = flipud(matrix_b); % Mirroring the matrix_b (flipping on horizontal axis) 
            Final_Data{k,5} = [matrix_c; matrix_a]; % Concatenate the two matrices and save at second position of the cell
       
        elseif K_locations_field{m}(n) > (Unprocessed_data(m).len - window) % Location K close to C terminus
            matrix_d = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)-window:Unprocessed_data(m).len,:); % Matrix containing pssm_prob from 32 downstream of K till end of protein seq
            matrix_e = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)-window:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Matrix containing pssm_prob downstream of K till the value to be mirrored
            matrix_f = flipud(matrix_e); % Mirroring the matrix_b (flipping on horizontal axis)
            Final_Data{k,5} = [matrix_d; matrix_f]; % Concatenate the two matrices and save at second position of the cell
            
        else % Location is good and does not require mirroring
            Final_Data{k,5} = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)-window:K_locations_field{m}(n)+window,:); % Save pssm_prob matrix (32 up and down streams of K) at second position
        
        end
        
        Final_Data{k,3} = Unprocessed_data(m).label{1}(K_locations_field{m}(n)); % Save class label at third position
        Final_Data{k,4} = K_locations_field{m}(n); % Save the K's location in the protein sequence

    end
   
end

Training_Data1 = Final_Data; 
clear Final_Data

% Bigram:

for l=1:size(Training_Data1,1)
    Training_Data1{l,2} = Training_Data1{l,5}/100; % Divide pssm_prob by 100
    
    Bigram_Mat = zeros(size(Training_Data1{l,5}, 2),size(Training_Data1{l,5}, 2));
    for M=1:size(Training_Data1{l,5}, 2)
       for N=1:size(Training_Data1{l,5}, 2)
           for i=1:size(Training_Data1{l,5},1)-1
            Bigram_Mat(M,N) = Bigram_Mat(M,N) + Training_Data1{l,5}(i,M)*Training_Data1{l,5}(i+1,N); % Bigram calculation
           end
       end
    end
    transposed_Bigram_Mat = Bigram_Mat'; % Transpose Bigram_Mat
    Training_Data1{l,2} = transposed_Bigram_Mat(:)'; % Store the transposed Bigram_Mat as column vector. This gives Feature Vector
end

Final_Data = Training_Data1;

%-------------------------------------------------------------
% Separation of neg and pos samples 
%-------------------------------------------------------------

Num_of_Pos  = 0;
Num_of_Neg  = 0;

for i = 1:size(Final_Data,1)
    
    if Final_Data{i,3} == '1'
       Num_of_Pos = Num_of_Pos + 1;
       Final_Pos(Num_of_Pos,:) = Final_Data(i,:); 
    end
    if Final_Data{i,3} == '0'
       Num_of_Neg = Num_of_Neg + 1;
       Neg_Samples(Num_of_Neg,:) = Final_Data(i,:);
    end    
end

%--------------------------------------------------------------------
% Filtering of neg samples to bring class imbalance from 1:29 to 1:2 
%--------------------------------------------------------------------

Filtered_Neg = Neg_Samples;

for a=1:size(Neg_Samples,1)
    
    for b=1:size(Final_Data, 1)

        Distance = pdist2(Neg_Samples{a,2}, Final_Data{b,2}); % Euclidean distance between the neg sample and a sample in whole data
   
        D(b, 1) = Distance; % Save the 3360 distances of entire data samples into 3360x1 matrix
        
    end
    
    for c=1:knn_neighbor % Number of nearest neighbors (including itself as it also calculated from itself at one point in time) to be taken into account for filtering. Number of neighbors is actually 111 
   
        [M,I] = min(D); % Find minimum from the 3360x1 matix and also its index
        
        if Final_Data{I,3} == '1' % If the minimum distance belongs to a +ve sample then delete the neg sample
            
            Filtered_Neg{a,2} = [];
            break
        end
       
        D(I) = [max(D)]; % Replace the minimum value with the maximum value of the matrix so that it wont point the same index as minimum in the next loop
        c = c+1;
    end
end

% Construct the neg samples after filtering
e = 1;
for j=1:size(Filtered_Neg,1)
    if isempty(Filtered_Neg{j,2}) == false
  
        Final_Neg(e,:) = Filtered_Neg(j,:); 
  
        e = e+1;
        
    end
end
 
%-----------------------------------------------------------------------------------
% Dividing data into 10 parts for 10 fold cross validation (these are 10 test sets)
%-----------------------------------------------------------------------------------

Filtered_Data_SamplesSep = [Final_Pos; Final_Neg]; % Combine the final sets of positive and negatives samples after the filtering process

len_of_data = size(Filtered_Data_SamplesSep, 1);

extra = mod(len_of_data,10);
base_size = (len_of_data - extra)/10;

for i = 1:10
    if (i>(10-extra))
        fold_size(i) = base_size + 1;
        
    else
        fold_size(i) = base_size;
        
    end  
end

default_gap = [5 10 15 20 25 30 35 40 45 50]; % Randomly assigned numbers (Number of +ve samples in each parts)
loop = 1;

% Dividing the Dataset into 10 parts and making sure +ve and -ve samples
% are well distributed
while(loop == 1)

        Data_Shuffled1 = Filtered_Data_SamplesSep(randperm(len_of_data),:);
        Data_Shuffled = Data_Shuffled1(randperm(len_of_data),:);

        Fold1_10 = Data_Shuffled(1:fold_size(1),:); % First Part
        j1=0;
        for i=1:fold_size(1)
        if strfind(Fold1_10{i,3},'1') == true
        j1 = j1+1;
        end
        end
        default_gap(1) = j1;

        Fold2_10 = Data_Shuffled(fold_size(1)+1:sum(fold_size(1:2)),:); % Second Part
        j2=0;
        for i=1:fold_size(2)
        if strfind(Fold2_10{i,3},'1') == true
        j2 = j2+1;
        end
        end
        default_gap(2) = j2;

        Fold3_10 = Data_Shuffled(sum(fold_size(1:2))+1:sum(fold_size(1:3)),:); % Third Part
        j3=0;
        for i=1:fold_size(3)
        if strfind(Fold3_10{i,3},'1') == true
        j3 = j3+1;
        end
        end
        default_gap(3) = j3;

        Fold4_10 = Data_Shuffled(sum(fold_size(1:3))+1:sum(fold_size(1:4)),:); % Fourth Part
        j4=0;
        for i=1:fold_size(4)
        if strfind(Fold4_10{i,3},'1') == true
        j4 = j4+1;
        end
        end
        default_gap(4) = j4;

        Fold5_10 = Data_Shuffled(sum(fold_size(1:4))+1:sum(fold_size(1:5)),:); % Fifth Part
        j5=0;
        for i=1:fold_size(5)
        if strfind(Fold5_10{i,3},'1') == true
        j5 = j5+1;
        end
        end
        default_gap(5) = j5;

        Fold6_10 = Data_Shuffled(sum(fold_size(1:5))+1:sum(fold_size(1:6)),:); % Sixth Part
        j6=0;
        for i=1:fold_size(6)
        if strfind(Fold6_10{i,3},'1') == true
        j6 = j6+1;
        end
        end
        default_gap(6) = j6;

        Fold7_10 = Data_Shuffled(sum(fold_size(1:6))+1:sum(fold_size(1:7)),:); % Seventh Part
        j7=0;
        for i=1:fold_size(7)
        if strfind(Fold7_10{i,3},'1') == true
        j7 = j7+1;
        end
        end
        default_gap(7) = j7;

        Fold8_10 = Data_Shuffled(sum(fold_size(1:7))+1:sum(fold_size(1:8)),:); % Eighth Part
        j8=0;
        for i=1:fold_size(8)
        if strfind(Fold8_10{i,3},'1') == true
        j8 = j8+1;
        end
        end
        default_gap(8) = j8;

        Fold9_10 = Data_Shuffled(sum(fold_size(1:8))+1:sum(fold_size(1:9)),:); % Ninth Part
        j9=0;
        for i=1:fold_size(9)
        if strfind(Fold9_10{i,3},'1') == true
        j9 = j9+1;
        end
        end
        default_gap(9) = j9;

        Fold10_10 = Data_Shuffled(sum(fold_size(1:9))+1:sum(fold_size(1:10)),:); % Tenth Part
        j10=0;
        for i=1:fold_size(10)
        if strfind(Fold10_10{i,3},'1') == true
        j10 = j10+1;
        end
        end
        default_gap(10) = j10;
        
        smallest = min(default_gap); % Minimum number of +ve samples in the overall parts
        largest = max(default_gap); % Maximum number of +ve samples in the overall parts
        
        if ((largest - smallest) < 3) % If difference between minimum and maximum number of +ve samples is less then 3 then there is good distribution of +ves and -ves
            
            loop = 0;
            
        end
end

%--------------------------------------------
% Contruction of training sets 
%--------------------------------------------

Whole_Data = {'Fold1_10' 'Fold2_10' 'Fold3_10' 'Fold4_10' 'Fold5_10' 'Fold6_10' 'Fold7_10' 'Fold8_10' 'Fold9_10' 'Fold10_10'}; 

for i = 1:10
     
    st = num2str(i);
    s = strcat('Fold', st, '_10');
    
    k = 1;
    for j = 1:10
        if regexp(Whole_Data{1,j}, s) == 1 % For each part, construct the corresponding unoptimized training set composed of other parts
            % Do nothing
        else
            index_train{1,k} = j;
            k = k+1;
        
        end
        
    end
    rt = num2str(index_train{1,1});
    r = strcat('Fold', rt, '_10');
    st = num2str(index_train{1,2});
    s = strcat('Fold', st, '_10');
    tt = num2str(index_train{1,3});
    t = strcat('Fold', tt, '_10');
    ut = num2str(index_train{1,4});
    u = strcat('Fold', ut, '_10');
    vt = num2str(index_train{1,5});
    v = strcat('Fold', vt, '_10');
    wt = num2str(index_train{1,6});
    w = strcat('Fold', wt, '_10');
    xt = num2str(index_train{1,7});
    x = strcat('Fold', xt, '_10');
    yt = num2str(index_train{1,8});
    y = strcat('Fold', yt, '_10');
    zt = num2str(index_train{1,9});
    z = strcat('Fold', zt, '_10');
    
    Train{i} = [eval(r); eval(s); eval(t); eval(u); eval(v); eval(w); eval(x); eval(y); eval(z)]; % Save the other 9 parts as unoptimized Training set
    
end

% Saving the training set into different variable

Train1 = Train{1};
Train2 = Train{2};
Train3 = Train{3};
Train4 = Train{4};
Train5 = Train{5};
Train6 = Train{6};
Train7 = Train{7};
Train8 = Train{8};
Train9 = Train{9};
Train10 = Train{10};
