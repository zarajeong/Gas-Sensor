%Optimasi Parameter

clear;
options = odeset('MaxStep',0.5,'RelTol',1e-6,'AbsTol',[1e-6,1e-6,1e-6,1e-6],'NonNegative',[1 2 3 4]);
    
%Suhu dalam Kelvin
Temp=165+273;
% Nilai awal (variasi)
nn=1;
Ef2=0;
Ef2part=0.005; 

% Looping 
while Ef2<=0.5
    %waktu 45 menit
         t=0;
         t0=t;
         tmax=45*60; %waktu akhir
    %inisialisasi variabel
         p=1;
         tpart=10; %loncatan
         TotalError=0; 
         TotalData=0;
         sigmaPt=1.31e+19;
    %Variabel coverage 
         x1=0;
         x2=0;
         x3=0;
         x4=0;
         
    %Input P diubah berdasarkan t
    while t<=tmax
          if and(t>=t0,t<t0+15*60) 
             PH2S=0/1E+6;
          elseif and(t>=t0+15*60,t<t0+30*60) 
                 PH2S=25/1E+6;
          elseif and(t>=t0+30*60,t<=t0+45*60) 
                 PH2S=0/1E+6;
          end

        %Nilai PO bergantung pada nilai H2S
        PO=0.20*(1-PH2S); %di data PO = 0 Pa
        T=Temp;
        
        % Penyelesaian PDB 
        [tt,x] = ode23tb(@(tt,x)PDB(tt,x,T,PH2S,PO,Ef2),[t t+tpart],[x1 x2 x3 x4],options);
        pmax = p+length(tt)-1;
        n=1;
      
        while p~=pmax 
              H2S(p)=x(n,1); 
              O(p)=x(n,2);
              SO2(p)=x(n,3); 
              H2O(p)=x(n,4);
              Total(p)=H2S(p)+O(p)+SO2(p)+H2O(p);
        
            % WF berbanding lurus dengan nilai coverage
            WFH2S(p)=H2S(p); 
            WFO(p)=O(p); 
            WFSO2(p)=SO2(p);
            WFH2O(p)=H2O(p);
            WFC(p)= -(WFH2S(p)+WFO(p)+WFSO2(p)+WFH2O(p))+0.0147;
       
              if and(tt(n)>=0,tt(n)<2710) 
                  Intercept=0.00144; 
                  a=1.0196E-4;
                  b=-7.97974E-7;
                  c=2.68259E-9;
                  d=-4.45246E-12;
                  e=4.01128E-15;
                  f=-2.07654E-18;
                  g=6.19422E-22;
                  h=-9.91165E-26;
                  i=6.59227E-30;
                  WFCeks(p)=Intercept + a*tt(n) + b*tt(n)^2 + c*tt(n)^3 + d*tt(n)^4 + e*tt(n)^5 + f*tt(n)^6 + g*tt(n)^7 + h*tt(n)^8 + i*tt(n)^9;
            end
            
           time(p)=tt(n);
           TotalError = TotalError+abs(WFC(p)-WFCeks(p));
           n=n+1;
           p=p+1;
        end

        p=p-1;
            t=t+tpart;
                t/45;
        x1=x(n-1,1);
        x2=x(n-1,2); 
        x3=x(n-1,3);
        x4=x(n-1,4);
    end

    AverageError(nn)=TotalError/p;
    figure; 
    EEf2(nn)=Ef2; 
        Ef2=Ef2+Ef2part;
             nn=nn+1;
% disp(TotalError)
% disp(AverageError(nn))
% disp(nn)
% disp(Ef1)
end

  plot(EEf2,AverageError,'c','LineWidth',1.25) %gcm
      xlabel('Ef2 (eV)');
          ylabel('Average Error (V)');
          
    % plot(time/60, WFCeks, 'g', time/60, WFC, 'r','LineWidth',1.25);
    %     xlabel('Time (m)');
    %        ylabel('CPD (V)');

 fid=fopen('Ef2x.txt', 'w'); 
      y=(EEf2); 
          fprintf(fid,'%20.18f\n',y);
              fclose(fid);
  fid=fopen('Ef2y.txt', 'w'); 
      y=(AverageError); 
          fprintf(fid,'%20.18f\n',y);
              fclose(fid);