%Dialog box with options in drop-down menu
function [choice,continueBool] = OptionsDialog(options,default_nr,title_str,prompt_str)
    
    %Open the dialog box
    Pix_SS = get(0,'screensize');
    middle = [Pix_SS(3)/2 Pix_SS(4)/2];
    boxSize = [Pix_SS(3)/5 Pix_SS(4)/5];
    d = dialog('Position',[middle(1)-Pix_SS(3)/10 middle(2)-Pix_SS(4)/10 boxSize],'Name',title_str);
    
    %Text
    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[boxSize(1)/7 boxSize(2)*5/7 boxSize(1)*5/7 boxSize(2)/7],...
           'String',prompt_str);
    
    %Drop-down menu   
    popup = uicontrol('Parent',d,...
           'Style','popup',...
           'Position',[boxSize(1)/7 boxSize(2)*3/7 boxSize(1)*5/7 boxSize(2)/7],...
           'String',options,...
           'Callback',@popup_callback);
    
    %Continue Button   
    btn = uicontrol('Parent',d,...
           'Position',[boxSize(1)/7 boxSize(2)*1/7 boxSize(1)*5/7 boxSize(2)/7],...
           'String','Continue',...
           'Callback',@btn_callback);
    
    %Set the default value   
    set(popup,'Value',default_nr);   
    choice = default_nr;
    continueBool = 0;
       
    %Wait for the dialog box to close before running to completion
    uiwait(d);
   
    %Callback function of the drop-down menu
    function popup_callback(popup,callbackdata)
       choice = get(popup,'Value');
    end

    %Callback function of the continue button
    function btn_callback(btn,callbackdata)
        continueBool = 1;
        delete(gcf);
    end
end