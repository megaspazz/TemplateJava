@ECHO OFF

SET SRC="src\#tools\generic.go"
SET EXE="bin\generic.exe"

SET FILE=%1.go

FOR /f "tokens=1,* delims= " %%a in ("%*") DO SET TYPES=%%b

go build -o %EXE% %SRC%
IF ERRORLEVEL 0 (
    %EXE% --file="%FILE%" --types="%TYPES%" --output="gen/_%FILE%"
)
