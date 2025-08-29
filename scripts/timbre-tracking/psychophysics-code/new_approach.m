function [COM] = new_approach(population_response_for_single_trial,CFs)

PROP = .9; %This is the proportion of the maximum peak height above which both peaks will considered of 'equal' height.
data = population_response_for_single_trial;
smoothed_data = supsmu(CFs,data,'Span',.5);         % Smooth data
smoothed_subtracted = max(data - smoothed_data,0);  % Subtract smoothed data


looking_for_maxes = sign(smoothed_subtracted);
contiguous_sections = runlength_function_Braden(looking_for_maxes);
section_no = 1;
section_indices = zeros(sum(contiguous_sections(1,:)),2);
for P = 1:size(contiguous_sections,2)
    if contiguous_sections(1,P) == 1
        if P ~= 1
            starting_index = sum(contiguous_sections(2,1:(P-1)))+1;
            ending_index = sum(contiguous_sections(2,1:P));
        elseif P == 1
            starting_index = 1;
            ending_index = contiguous_sections(2,1);
        end
        section_indices(section_no,1:2) = [starting_index,ending_index];
        section_no = section_no + 1;
    end
end
maxes = [];
for Q = 1:(section_no-1)
    maxes(Q) = max(data(section_indices(Q,1):section_indices(Q,2)));
end
sortable_maxes = [maxes;1:length(maxes)]';
sortable_maxes2 = sortrows(sortable_maxes,'descend');
two_sections = sortable_maxes2(1:2,2);
indices_section_1 = section_indices(two_sections(1),:);
indices_section_2 = section_indices(two_sections(2),:);
data_section_1 = data(indices_section_1(1):indices_section_1(2));
data_section_2 = data(indices_section_2(1):indices_section_2(2));

one_split_peak = sign( mean(data_section_2)./mean(data_section_1) - PROP); 
final_vector_3 = zeros(1,length(CFs));
if one_split_peak == 1 || one_split_peak == 0
    final_vector_3(indices_section_1(1):indices_section_1(2)) = data_section_1; %Using data here instead of slope-subtracted to involve local thresholding less.
    final_vector_3(indices_section_2(1):indices_section_2(2)) = data_section_2; %Using data here instead of slope-subtracted to involve local thresholding less.
elseif one_split_peak == -1
    final_vector_3(indices_section_1(1):indices_section_1(2)) = data_section_1; %Using data here instead of slope-subtracted to involve local thresholding less.
end
final_vector = final_vector_3;

COM = sum(CFs.*final_vector)/sum(final_vector);
end


function [encoded] = runlength_function_Braden(looking_for_maxes)

for OO = 1:length(looking_for_maxes)
    if OO == 1
        values = looking_for_maxes(OO);
        runlengths = 1;
    else
        if looking_for_maxes(OO) == looking_for_maxes(OO-1)
            runlengths(end) = runlengths(end)+1;
        elseif looking_for_maxes(OO) ~= looking_for_maxes(OO-1)
            values = [values looking_for_maxes(OO)];
            runlengths = [runlengths 1];
        end
    end
end

encoded = [values; runlengths];
end