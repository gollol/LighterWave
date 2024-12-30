@echo off

SET SCRIPT_DIR=%~dp0
SET SRC_DIR=%~dp0..\

SET BOB_DIR=%SRC_DIR%\deps\bob
SET LIB_DIR=%SCRIPT_DIR%\lightwave_blender

cd %BOB_DIR%
py -3.11 setup.py build_ext --build-lib=%LIB_DIR%

cd %SCRIPT_DIR%
tar.exe -a -cf lightwave_blender.zip lightwave_blender
if %errorlevel% neq 0 exit /b %errorlevel%


echo Blender plugin built. Open Blender and go to 'Edit - Preferences - Addons - Install...'

echo Browse to the 'lightwave_blender.zip' file in this directory and install it.