function getParameters(classifier,cv_setting)
%
% This function retrieves optimal parameter values based on the supplied
% algorithm, the CV scenario used and the dataset on which prediction will
% be performed.
%
% INPUT:
%  classifier:  algorithm to be used
%  cv_setting:          
%
% OUTPUT:
%  params:   optimal parameter values for given inputs


global k alpha  k1 k2  pp lamda mu1 mu2    
global num_iter p lambda_l lambda_d lambda_t    


    
    switch classifier
      
      case 'mf'
                    switch cv_setting 
                        case 1
                            k=5;alpha=1.5;    
                          case 2
                            k=15;alpha=0.01; 
                        case 3
                            k=20;alpha=0.01; 
                      end
       
      case 'dmf'
          switch cv_setting 
                        case 1
                            alpha=2; k1=25;k2=5;    
                          case 2
                            alpha=0.5; k1=25;k2=20;
                          case 3
                            alpha=0.5; k1=30;k2=15;
                      end
            
            
      case 'grmf'
          
            num_iter = 2;
            k = 100;
                      switch(cv_setting)
                        
                        case 1
                           p=7; lambda_l = 0.25; lambda_d = 0.25;lambda_t = 0.05;    
                          case 2
                           p=3; lambda_l = 0.125; lambda_d = 0.3;lambda_t = 0.1; 
                           case 3
                           p=; lambda_l = ; lambda_d = ;lambda_t = ; 
                      end
                                   
                                   
        case 'gr1bmc_ppxa'
           
                    switch cv_setting 
                        
                        case 1
                            pp=2;lamda=0.1;mu1=0.1;mu2=0.1;
                        case 2
                            pp=2;lamda=0.01;mu1=1;mu2=1;
                       case 3
                            %pp=;lamda=;mu1=;mu2=;
                    end
        case 'dgrdmf'
            switch cv_setting 
                        
                case 1
                     pp=8;alpha=2;k1=25;k2=15;mu1=100;mu2=0.1;
                case 2
                    pp=2;alpha=2;k1=35;k2=10;mu1=100;mu2=0.01;
                case 3
                    %pp=;alpha=;k1=;k2=;mu1=;mu2=;
            end
            
        case 'grmc_admm'
                    switch cv_setting 
                        case 1
                            pp=2;lamda=0.01;mu1=0.01;mu2=0.005;
                        case 2
                            pp=2;lamda=0.01;mu1=0.1;mu2=0.01;
                        case 3
                            %pp=;lamda=;mu1=;mu2=;
                    end
        case 'mc'
              
    end
            

end