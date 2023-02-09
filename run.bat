@ECHO OFF

SETLOCAL EnableDelayedExpansion

SET fullpath=%~1
SET filename=%~n1
CALL SET filepath=%%fullpath:%~dp0=%%

SET JAVA_BIN=C:\Program Files\Java\jdk1.8.0_361\bin
SET JAVA=java.exe
SET JAVAC=javac.exe

IF NOT [%JAVA_BIN] == [] (
    SET JAVA=%JAVA_BIN%\%JAVA%
    SET JAVAC=%JAVA_BIN%\%JAVAC%
)

SET COMPILE="%JAVAC%" -Xlint -d "bin" "%filepath%"
SET EXECUTE="%JAVA%" -Xmx1024m -Xss256m -cp bin "%filename%"

SET INPUT=io\in.txt
SET OUTPUT=io\out.txt
SET ERROR=io\err.txt

IF NOT EXIST "bin\" (
    MKDIR "bin"
)

ECHO.
IF DEFINED filepath (
    ECHO === Compiling:  %filepath%
    %COMPILE% 2> %ERROR%
    IF ERRORLEVEL 1 (
        ECHO.
        TYPE "%ERROR%"
        ECHO.
        ECHO === Compilation failed.  See "%ERROR%" for details.
    ) ELSE (
        FOR %%A IN (%ERROR%) DO IF NOT %%~zA == 0 SET hasCompilationErrorOutput=true
        ECHO.
        IF DEFINED hasCompilationErrorOutput (
            TYPE "%ERROR%"
            ECHO.
            ECHO === Compilation succeeded with warnings.
        ) ELSE (
            ECHO === Compilation successful.
        )
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
        )
    )
) ELSE (
    ECHO === ERROR: File not found!
)
ECHO.

ENDLOCAL