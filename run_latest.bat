@ECHO OFF

SETLOCAL EnableDelayedExpansion

SET SEARCH_ROOT=src
SET FILE_EXTENSION=java

FOR /f "usebackq delims=" %%a IN (`
  POWERSHELL -NoP -C "(Get-ChildItem -Path '%SEARCH_ROOT%' -Filter '*.%FILE_EXTENSION%' -File -Recurse | Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName"
`) DO SET "filepath=%%a"

CALL .\run.bat "%filepath%"

ENDLOCAL