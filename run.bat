@ECHO OFF

SETLOCAL EnableDelayedExpansion

FOR /f "delims=" %%F IN ('dir /b /s "src\%1.java" 2^>NUL') DO SET filepath=%%F
FOR %%F IN (%filepath%) DO SET filename=%%~nF

SET INPUT=io\in.txt
SET OUTPUT=io\out.txt
SET ERROR=io\err.txt

SET COMPILE=javac "%filepath%" -d "bin"
SET EXECUTE=java -cp bin "%filename%"

IF NOT EXIST "bin\" (
    MKDIR "bin"
)

IF DEFINED filepath (
    ECHO === Compiling:  %filepath%
    %COMPILE% 2> %ERROR%
    IF ERRORLEVEL 1 (
        ECHO === Compilation failed.  See "%ERROR%" for details.
        ECHO.
        TYPE "%ERROR%"
        ECHO.
    ) ELSE (
        ECHO === Compilation successful.
        ECHO.
        %EXECUTE% < %INPUT% > %OUTPUT% 2> %ERROR% ^
            && SET success=true ^
            || SET success=
        TYPE "%OUTPUT%"
        ECHO.
        if DEFINED success (
            ECHO === Execution successful.
        ) ELSE (
            ECHO === The program crashed. See "%ERROR%" for details.
            ECHO.
            TYPE "%ERROR%"
            ECHO.
        )
    )
) ELSE (
    ECHO === ERROR: File not found!
)

ENDLOCAL
