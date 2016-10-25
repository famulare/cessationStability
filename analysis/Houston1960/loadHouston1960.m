function [index,sibling,contact]=loadHouston1960()

for k=1:3
    index(k)=xls2struct('..\..\Data\Houston1960\index.xlsx',['mOPV',num2str(k)]);
    sibling(k)=xls2struct('..\..\Data\Houston1960\sibling.xlsx',['mOPV',num2str(k)]);
    contact(k)=xls2struct('..\..\Data\Houston1960\primaryContact.xlsx',['mOPV',num2str(k)]);
end

for k=1:3
    % contacts
        % age correction factor for contacts
        % contacts have similar demographics to siblings, but demographic
        % breakdown is not presented in Benyesh-Melnick1967 paper
        estimatedDenominatorFractionUnder5 = (sum([sibling.DenominatorUnder5])./sum([sibling.DenominatorTotal]));
        estimatedSheddingFractionUnder5 = (sum([sibling.SheddingUnder5])./sum([sibling.SheddingTotal]));
        
        contact(k).DenominatorUnder5=round(contact(k).DenominatorTotal*estimatedDenominatorFractionUnder5(k));
        contact(k).SheddingUnder5=round(contact(k).SheddingTotal*estimatedSheddingFractionUnder5(k));
        contact(k).Denominator5plus=round(contact(k).DenominatorTotal*(1-estimatedDenominatorFractionUnder5(k)));
        contact(k).Shedding5plus=round(contact(k).SheddingTotal*(1-estimatedSheddingFractionUnder5(k)));
end

end
