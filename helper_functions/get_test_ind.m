function [test_ind,len] = get_test_ind( cv_setting,i,n, Y,Sd,Sv)
%GET_TEST_IND Summary of this function goes here
%   Detailed explanation goes here

rng(0);
num_drugs=size(Sd,1); num_virus=size(Sv,1);

       if(cv_setting==1)
       %CV1
           len=numel(Y); 
           rand_ind = randperm(len); 
           test_ind = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
       
           
       elseif(cv_setting==2)
       %CV2     
           
           len=num_virus; rand_ind = (randperm(len));
           left_out_virus = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
           test_ind = zeros(num_drugs,length(left_out_virus));
           for j=1:length(left_out_virus)
                    curr_left_out_virus = left_out_virus(j);
                    %test_ind(j,:) = ((0:(num_drugs-1)) .* num_virus) + curr_left_out_virus;
                    test_ind(:,j) = (1:num_drugs)' + ((curr_left_out_virus-1)*num_drugs);
           end
           test_ind = reshape(test_ind,numel(test_ind),1);
            
           
       elseif(cv_setting==3)
           len=num_drugs; rand_ind = (randperm(len));
           left_out_drugs = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
            test_ind = zeros(length(left_out_drugs),num_virus);
            for j=1:length(left_out_drugs)
                curr_left_out_drug = left_out_drugs(j);
                test_ind(j,:) = ((0:(num_virus-1)) .* num_drugs) + curr_left_out_drug;
            end
            test_ind = reshape(test_ind,numel(test_ind),1);
            
       end
       test_ind=test_ind(:);
end



