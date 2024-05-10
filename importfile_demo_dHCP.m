function demodHCPrev = importfile(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  DEMODHCPREV = IMPORTFILE(FILE) reads data from the first worksheet in
%  the Microsoft Excel spreadsheet file named FILE.  Returns the data as
%  a table.
%
%  DEMODHCPREV = IMPORTFILE(FILE, SHEET) reads from the specified
%  worksheet.
%
%  DEMODHCPREV = IMPORTFILE(FILE, SHEET, DATALINES) reads from the
%  specified worksheet for the specified row interval(s). Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  demodHCPrev = importfile("/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/05_analysis_cmap_pmap/000_sourceData_Code/Upload_versions/files/demo_dHCP_rev.xlsx", "Sheet1", [2, 375]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 2024-04-16 11:10:56

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 375];
end

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":L" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["sort", "dataset", "sub_id", "ses_id", "scanage", "birthage", "sex", "radscore", "preproc", "meanFD", "MeanFD04", "qc_grad"];
opts.VariableTypes = ["double", "categorical", "string", "double", "double", "double", "categorical", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "sub_id", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["dataset", "sub_id", "sex"], "EmptyFieldRule", "auto");

% Import the data
demodHCPrev = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":L" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    demodHCPrev = [demodHCPrev; tb]; %#ok<AGROW>
end

end