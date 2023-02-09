@ECHO OFF

SETLOCAL EnableDelayedExpansion

SET SEARCH_ROOT="src"
SET FILE_EXTENSION=java

SET SEARCH_PATH="%SEARCH_ROOT%\%1.%FILE_EXTENSION%"

FOR /f "delims=" %%F IN ('dir /b /s "%SEARCH_PATH%" 2^>NUL') DO SET filepath=%%F

CALL .\run.bat "%filepath%"

ENDLOCAL