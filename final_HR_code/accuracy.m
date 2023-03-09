function [mean_accuracy, median_accuracy, accuracy_values]=accuracy(evolutionHR,evolutionReference)
% calculate accuracy

if length(evolutionHR) == length(evolutionReference)
    accuracy_values = zeros(length(evolutionHR),1);
    for i =1:length(evolutionHR)
        accuracy_values(i)=((evolutionReference(i)-abs(evolutionReference(i)-evolutionHR(i)))/evolutionReference(i))*100;
    end

    mean_accuracy = mean(accuracy_values);
    median_accuracy = median(accuracy_values);

else
    mean_accuracy = 0;
    median_accuracy = 0;
end
end

