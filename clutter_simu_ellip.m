%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clutter_out=clutter_simu_ellip(area_num,simu_range,simu_num,data_len,clut_num,para)
%assume the spatial distribution of each high clutter density area is Gauss
%area_num, the number of high clutter area
%simu_range, the total area of surveillance area, assume its shape is rectangular, x_min, x_max,y_min,y_max
%simu_num, number of simulation experiments
%data_len, total sample in one smiulation experiments
%clut_num, expected clutter number in whole surveillance region for each scan
%par, parameter for each high clutter area, its format :
%(1) each row corresponding to one high clutter area;
%(2) for each row: percentenge of clutter number,x_center,x_sigma,y_center,y_sigma
%(3) row number of para should be equal to area_num
%clutter_out, 2-D cell matrix, each cell is a 2 row matrix, which has contained clutter's position in one sample frame 
clutter_out=cell(simu_num,data_len);
for k=1:simu_num
    for i=1:data_len
        temp_data=zeros(1,1);
        dump_index=1;
        all_num=random('poiss',clut_num);
        all_num=2*clut_num;
        while (all_num>=2*clut_num)|(all_num<=0.5*clut_num)
            all_num=random('poiss',clut_num);
        end
        if area_num>=1
            for j=1:area_num
                temp_num=round(all_num*para(j,1));
                if temp_num>=1
                    for m=1:temp_num
                        accetp_flag=0;
                        while accetp_flag==0
                            dump=para(j,2)+para(j,3)*randn;
                            if (dump>=simu_range(1))&(dump<=simu_range(2))
                                temp_data(:,dump_index)=dump;
                                dump_index=dump_index+1;
                                accetp_flag=1;
                            end
                        end
                    end
                end
            end
            total_high=sum(para(:,1));
            if total_high<1
                temp_num=round(all_num*(1-total_high));
                if temp_num>=1
                    for m=1:temp_num
                        temp_data(:,dump_index)=simu_range(1)+(simu_range(2)-simu_range(1))*rand;
                        dump_index=dump_index+1;
                    end
                end
            end
        else
            temp_num=random('poiss',clut_num);
            if temp_num>=1
                for m=1:temp_num
                    temp_data(:,dump_index)=simu_range(1)+(simu_range(2)-simu_range(1))*rand;
                    dump_index=dump_index+1;
                end
            end
        end
        clutter_out{k,i}=temp_data;
    end
end