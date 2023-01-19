classdef db_eval_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        SnapshotButton           matlab.ui.control.Button
        ZumEditField             matlab.ui.control.NumericEditField
        ZumEditFieldLabel        matlab.ui.control.Label
        ReadDBFileButton         matlab.ui.control.Button
        PropertiesTextArea       matlab.ui.control.TextArea
        PropertiesTextAreaLabel  matlab.ui.control.Label
        ZumSlider                matlab.ui.control.Slider
        ZumSliderLabel           matlab.ui.control.Label
        UIAxes6                  matlab.ui.control.UIAxes
        UIAxes5                  matlab.ui.control.UIAxes
        UIAxes4                  matlab.ui.control.UIAxes
        UIAxes3                  matlab.ui.control.UIAxes
        UIAxes2                  matlab.ui.control.UIAxes
        UIAxes                   matlab.ui.control.UIAxes
    end


    properties (Access = private)
        z; % Description
        effic;
        efficsh;
        fileSTR;
        t;
        nu;
        c0 = 3e8;
        lambda;
        fileDirs;

        %% Plots
        plot1;
        plot2;
        plot3;
        plot4;
        plot5;
        plot6;
        plotEfficMax;
        plotEfficCurrent;
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ReadDBFileButton
        function ReadDBFileButtonPushed(app, event)
            [fileName,fileLoc] = uigetfile("*.*","MultiSelect","off");
            if fileName ==0; return; end
            clc;
            %c0 = 3e8;
            app.fileSTR = strcat(fileLoc,fileName);
            app.fileDirs = fileLoc;
            app.z = h5read(app.fileSTR,"/z")*1e6;
            app.t = h5read(app.fileSTR,"/t")*1e12;
            app.nu= h5read(app.fileSTR,"/nu")*1e-12;
            app.lambda = [0,app.c0./app.nu(2:end)]*1e-6;
            app.ZumSlider.Limits = [min(app.z),max(app.z)];
            app.effic = h5read(app.fileSTR,"/effic")*100;
            app.plot1 = plot(app.UIAxes, app.z,app.effic);
            [M,I] = max(app.effic);
            ETHz = h5read(app.fileSTR,"/"+num2str(app.z(I))+"/ETHz")/1e5;
            ATHz = normalize(h5read(app.fileSTR,"/"+num2str(app.z(I))+"/ATHz"),"range");
            Eop = h5read(app.fileSTR,"/"+num2str(app.z(I))+"/Eop")/1e13;
            Aop = normalize(h5read(app.fileSTR,"/"+num2str(app.z(I))+"/Aop"),"range");
            app.plot2 = plot(app.UIAxes2,app.t,ETHz);
            xlim(app.UIAxes2,[-3,3]);
            app.plot3 = plot(app.UIAxes3,app.nu,ATHz);
            app.plot5 = plot(app.UIAxes5,app.t,Eop);
            app.plot6 = plot(app.UIAxes6,app.lambda,Aop);
            lambdaCentral = app.c0/sum(Aop(2:end).*app.nu(2:end).')*sum(Aop(2:end))*1e-6;
            xlim(app.UIAxes3,[0,5]);
            xlim(app.UIAxes5,[-3,3]);
            xlim(app.UIAxes6,[lambdaCentral-3,lambdaCentral+3]);
            app.ZumSlider.Value = app.z(I);
            app.PropertiesTextArea.Value{1} = ['Maximum efficiency = ',num2str(M),' @ ', num2str(app.z(I)), 'um'];
            app.PropertiesTextArea.Value{2} = ['Maximum THz-field = ',num2str(max(abs(ETHz))),' kV/cm', ' @ ', num2str(app.z(I)), 'um'];
            app.PropertiesTextArea.Value{3} = ['Efficiency @ ', num2str(app.z(I)), 'um = ', num2str(app.effic(I)), '%'];
            app.PropertiesTextArea.Value{4} = ['THz central frequency @ ',num2str(app.z(I)), ' um = ' num2str(sum(ATHz.*app.nu.')/sum(ATHz)), ' THz'];
            app.PropertiesTextArea.Value{5} = ['Maximum pump intensity @ ', num2str(app.z(I)),' = ', num2str(max(Eop)),' GW/cm^2'];
            app.PropertiesTextArea.Value{6} = ['Pump central wavelength @ ',num2str(app.z(I)), ' um = ' num2str(lambdaCentral), ' um'];
            hold(app.UIAxes,"on")
            app.plotEfficMax = plot(app.UIAxes,app.z(I),app.effic(I),"r*");
            app.plotEfficCurrent = plot(app.UIAxes,app.z(I),app.effic(I),"g*");
            hold(app.UIAxes,"off")
        end

        % Callback function: ZumSlider, ZumSlider
        function ZumSliderValueChanged(app, event)
            value = event.Value;
            [~,z_index] = min(abs(app.z-value));
            app.ZumSlider.Value=app.z(z_index);
            ETHz = h5read(app.fileSTR,"/"+num2str(app.z(z_index))+"/ETHz")/1e5;
            ATHz = normalize(h5read(app.fileSTR,"/"+num2str(app.z(z_index))+"/ATHz"),"range");
            Eop = h5read(app.fileSTR,"/"+num2str(app.z(z_index))+"/Eop")/1e13;
            Aop = normalize(h5read(app.fileSTR,"/"+num2str(app.z(z_index))+"/Aop"),"range");
            app.plot2.YData = ETHz;
            app.plot3.YData = ATHz;
            app.plot5.YData = Eop;
            app.plot6.YData = Aop;
            lambdaCentral = app.c0/sum(Aop(2:end).*app.nu(2:end).')*sum(Aop(2:end))*1e-6;
            xlim(app.UIAxes6,[lambdaCentral-3,lambdaCentral+3]);
            app.plotEfficCurrent.XData = app.z(z_index);
            app.plotEfficCurrent.YData = app.effic(z_index);
            app.PropertiesTextArea.Value{2}= ['Maximum THz-field = ',num2str(max(abs(ETHz))),' kV/cm', ' @ ', num2str(app.z(z_index)), 'um'];
            app.PropertiesTextArea.Value{3} = ['Efficiency @ ', num2str(app.z(z_index)), 'um = ', num2str(app.effic(z_index)), '%'];
            app.PropertiesTextArea.Value{4} = ['THz central frequency @ ',num2str(app.z(z_index)), ' um = ' num2str(sum(ATHz.*app.nu.')/sum(ATHz)), 'THz'];
            app.PropertiesTextArea.Value{5} = ['Maximum pump intensity @ ', num2str(app.z(z_index)),' = ', num2str(max(Eop)),' GW/cm^2'];
            app.PropertiesTextArea.Value{6} = ['Pump central wavelength @ ',num2str(app.z(z_index)), ' um = ' num2str(lambdaCentral), ' um'];
            app.ZumEditField.Value = app.ZumSlider.Value;
        end

        % Value changed function: ZumEditField
        function ZumEditFieldValueChanged(app, event)
            value = app.ZumEditField.Value;
            app.ZumSlider.Value = value;
            ZumSliderValueChanged(app,event)
            app.ZumEditField.Value = app.ZumSlider.Value;
        end

        % Button pushed function: SnapshotButton
        function SnapshotButtonPushed(app, event)
            z_current = app.ZumSlider.Value;
            ETHz = h5read(app.fileSTR,"/"+num2str(z_current)+"/ETHz")/1e5;
            ATHz = normalize(h5read(app.fileSTR,"/"+num2str(z_current)+"/ATHz"),"range");
            Eop = h5read(app.fileSTR,"/"+num2str(z_current)+"/Eop")/1e13;
            Aop = normalize(h5read(app.fileSTR,"/"+num2str(z_current)+"/Aop"),"range");
            writematrix([app.z(:),app.effic(:)],[app.fileDirs,'/effic.txt']);
            dir = [app.fileDirs,'/Snapshot_',num2str(z_current),'_um'];
            mkdir(dir);
            writematrix([app.t(:),ETHz(:)],[dir,'/ETHz']);
            writematrix([app.t(:),Eop(:)],[dir,'/Eop']);
            writematrix([app.nu(:),ATHz(:)],[dir,'/ATHz']);
            writematrix([app.lambda(:),Aop(:)],[dir,'/Aop']);
            writecell(app.PropertiesTextArea.Value,[dir,'/props']);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1003 561];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Conversion Efficiency')
            xlabel(app.UIAxes, 'Crystal Length (mm)')
            ylabel(app.UIAxes, 'Conversion Efficiency (%)')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Box = 'on';
            app.UIAxes.Position = [54 363 300 185];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'THz Field')
            xlabel(app.UIAxes2, 'Time (ps)')
            ylabel(app.UIAxes2, 'Field Strength (kV/cm)')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Box = 'on';
            app.UIAxes2.Position = [354 363 300 185];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.UIFigure);
            title(app.UIAxes3, 'THz Spectrum')
            xlabel(app.UIAxes3, 'Frequency (THz)')
            ylabel(app.UIAxes3, 'Spectral Amplitude')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Box = 'on';
            app.UIAxes3.Position = [654 363 300 185];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.UIFigure);
            title(app.UIAxes4, 'THz Spectrum')
            xlabel(app.UIAxes4, 'Crystal Length (mm)')
            ylabel(app.UIAxes4, 'Conversion efficiency (%)')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.Box = 'on';
            app.UIAxes4.Position = [55 179 300 185];

            % Create UIAxes5
            app.UIAxes5 = uiaxes(app.UIFigure);
            title(app.UIAxes5, 'Pump Intensity')
            xlabel(app.UIAxes5, 'Time (ps)')
            ylabel(app.UIAxes5, 'Intensity (GW/cm^2)')
            zlabel(app.UIAxes5, 'Z')
            app.UIAxes5.Box = 'on';
            app.UIAxes5.Position = [355 179 300 185];

            % Create UIAxes6
            app.UIAxes6 = uiaxes(app.UIFigure);
            title(app.UIAxes6, 'Pump Spectrum')
            xlabel(app.UIAxes6, 'Wavelength (um)')
            ylabel(app.UIAxes6, 'Spectral Amplitude')
            zlabel(app.UIAxes6, 'Z')
            app.UIAxes6.Box = 'on';
            app.UIAxes6.Position = [654 179 300 185];

            % Create ZumSliderLabel
            app.ZumSliderLabel = uilabel(app.UIFigure);
            app.ZumSliderLabel.HorizontalAlignment = 'right';
            app.ZumSliderLabel.Position = [509 125 40 22];
            app.ZumSliderLabel.Text = 'Z (um)';

            % Create ZumSlider
            app.ZumSlider = uislider(app.UIFigure);
            app.ZumSlider.ValueChangedFcn = createCallbackFcn(app, @ZumSliderValueChanged, true);
            app.ZumSlider.ValueChangingFcn = createCallbackFcn(app, @ZumSliderValueChanged, true);
            app.ZumSlider.Position = [570 134 374 3];

            % Create PropertiesTextAreaLabel
            app.PropertiesTextAreaLabel = uilabel(app.UIFigure);
            app.PropertiesTextAreaLabel.HorizontalAlignment = 'right';
            app.PropertiesTextAreaLabel.Position = [38 118 60 22];
            app.PropertiesTextAreaLabel.Text = 'Properties';

            % Create PropertiesTextArea
            app.PropertiesTextArea = uitextarea(app.UIFigure);
            app.PropertiesTextArea.Position = [113 48 325 94];

            % Create ReadDBFileButton
            app.ReadDBFileButton = uibutton(app.UIFigure, 'push');
            app.ReadDBFileButton.ButtonPushedFcn = createCallbackFcn(app, @ReadDBFileButtonPushed, true);
            app.ReadDBFileButton.Position = [855 36 100 23];
            app.ReadDBFileButton.Text = 'Read DB File';

            % Create ZumEditFieldLabel
            app.ZumEditFieldLabel = uilabel(app.UIFigure);
            app.ZumEditFieldLabel.HorizontalAlignment = 'right';
            app.ZumEditFieldLabel.Position = [678 36 44 22];
            app.ZumEditFieldLabel.Text = 'Z  (um)';

            % Create ZumEditField
            app.ZumEditField = uieditfield(app.UIFigure, 'numeric');
            app.ZumEditField.ValueChangedFcn = createCallbackFcn(app, @ZumEditFieldValueChanged, true);
            app.ZumEditField.Position = [737 36 100 22];

            % Create SnapshotButton
            app.SnapshotButton = uibutton(app.UIFigure, 'push');
            app.SnapshotButton.ButtonPushedFcn = createCallbackFcn(app, @SnapshotButtonPushed, true);
            app.SnapshotButton.Position = [559 35 100 23];
            app.SnapshotButton.Text = 'Snapshot';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = db_eval_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end