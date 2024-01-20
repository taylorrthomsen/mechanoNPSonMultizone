% function for choosing thresholds

% if isempty(thresholds)
%     disp("please enter thresholds value for ref and queeze into command window")
% end

%inputs = 
%thresholds vector

%outputs = 
%thresholds vector

function [thresholds] = ...
    thresholdingTT(diff1,thresholds, cellCount, startidx)


keepGoing = 1;
%userinput = 1
% is user doesnt want to change the thing, just use enter and the thing
% will become empty
    
% figure();
    
while keepGoing == 1
%      fa = gcf;
%     ax = gca;
    gca;
%     hold on
    plot(diff1, 'color', 'b'); 
    hold on; yline(-thresholds(1),'color','r'); yline(-thresholds(2),'color','g'); yline(thresholds(1),'color','r'); yline(thresholds(2),'color','g')
%     titleString = sprintf('%s - %s', cellCount, startidx);
    title(cellCount, startidx)
    hold off
    
%     userInputSkip = input("Would you like to skip this event? (1 for Y/2 for N)");
% 
%      if userInputSkip == 1
%         break;
%      end

    userInputRef = input("Would you like to change the ref threshold? (1 for Y/2 for N)");
    userInputSq = input("Would you like to change the squeeze threshold? (1 for Y/2 for N)");


    if (userInputRef == 1) || (userInputSq == 1)
     if userInputRef == 1
        newRef = input("type new ref val");
        thresholds(1) = newRef;
     end
     if userInputSq == 1
        newSq = input("type new squeeze val");
        thresholds(2) = newSq;
     end

    else 
        keepGoing=2;
    end

end
end



%     figure;
%     
%     plot(diff1, 'color', 'b'); hold on; yline(-thresholds(1),'color','r'); hold on; yline(-thresholds(2),'color','g'); hold on; yline(thresholds(1),'color','r'); hold on; yline(thresholds(2),'color','g')
%     hold off

    %for after

%    userInput1 = input('Please input new ref threshold If you want to keep current ref thresold, enter 0:');
%    if ((userInput1 ~= thresholds(2)) && (userInput1 ~= 0))
%        thresholds(1) = userInput1;
%    elseif userInput1 == 0
%         
%    else
%         fprintf("please enter new threhold or enter 0")
%    end
% 
%    userInput2 = input('Please input new squeeze threshold If you want to keep current squeeze thresold, enter 0:');
%    if ((userInput2 ~= thresholds(2)) && (userInput2 ~= 0))
%        thresholds(2) = userInput2;
%    elseif userInput2 == 0
%         
%    else
%         fprintf("please enter new threhold or enter 0")
%    end
% 
%    figure;
%    plot(diff1, 'color', 'b'); hold on; yline(-thresholds(1),'color','r'); hold on; yline(-thresholds(2),'color','g'); hold on; yline(thresholds(1),'color','r'); hold on; yline(thresholds(2),'color','g')
%    hold off
% 
%    userInput3 = input('Please input new ref threshold If you want to keep current ref thresold, enter 0:');
%    if ((userInput3 ~= thresholds(2)) && (userInput3 ~= 0))
%        thresholds(1) = userInput3;
%    elseif userInput3 == 0
%         
%    else
%         fprintf("please enter new threhold or enter 0")
%    end
%  
%     userInput4 = input('Please input new squeeze threshold If you want to keep current squeeze thresold, enter 0:');
%     if ((userInput4 ~= thresholds(2)) && (userInput4 ~= 0))
%         thresholds(2) = userInput4;
%     elseif userInput4 == 0
%          
%     else
%          fprintf("please enter new threhold or enter 0")
%     end
%     
%    figure;
%    plot(diff1, 'color', 'b'); hold on; yline(-thresholds(1),'color','r'); hold on; yline(-thresholds(2),'color','g'); hold on; yline(thresholds(1),'color','r'); hold on; yline(thresholds(2),'color','g')
%    hold off

  % eventEnd =yas2det((locs(i))-5300);
  % find(yas2det,eventStart)

