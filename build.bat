
@echo off

set INTERMEDIATE=.\intermediate
set DEBUG_RUN_TREE=.\run_trees\debug
set EXE=engine.exe
set PDB=engine.pdb

set SOURCE=src\unity_build.cpp
REM set LIBRARY_SOURCE=lib\glad\src\glad.c lib\imgui\imgui*.cpp lib\imgui\examples\imgui_impl_glfw.cpp lib\imgui\examples\imgui_impl_opengl3.cpp
set LIBRARY_SOURCE=
set INCLUDE_DIRS=/I"src" /I"lib\glad\include" /I"lib\glfw\include" /I"lib\imgui" /I"lib\imgui\examples" /I"lib\stb"
set LIBS=user32.lib gdi32.lib shell32.lib opengl32.lib lib\glfw\glfw.lib

set COMMON_COMPILE_FLAGS=/c /EHsc /std:c++20 /Fo%INTERMEDIATE%\

set DEBUG_MACROS=/DDEBUG
set DEBUG_COMPILE_FLAGS=/Zi
set DEBUG_LINK_FLAGS=/DEBUG:FULL /OUT:"%INTERMEDIATE%\%EXE%"

mkdir %INTERMEDIATE%
mkdir %DEBUG_RUN_TREE%

cl %COMMON_COMPILE_FLAGS% %DEBUG_COMPILE_FLAGS% %DEBUG_MACROS% %INCLUDE_DIRS% %SOURCE% %LIBRARY_SOURCE%
if %errorlevel% neq 0 exit /b %errorlevel%
link %DEBUG_LINK_FLAGS% %LIBS% %INTERMEDIATE%\*.obj
if %errorlevel% neq 0 exit /b %errorlevel%

xcopy /Y %INTERMEDIATE%\%EXE% %DEBUG_RUN_TREE%
xcopy /Y %INTERMEDIATE%\%PDB% %DEBUG_RUN_TREE%
xcopy /Y /E assets %DEBUG_RUN_TREE%\assets\



echo.
echo %DEBUG_RUN_TREE%\%EXE%
echo.
REM %DEBUG_RUN_TREE%\%EXE%

