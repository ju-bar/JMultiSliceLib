@echo off
echo ---------------------------------------------------
echo Updating JMultiSliceLib C++ Binary Backups.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
echo Removing previous zip files.
del "JMSLib_BIN_Backup.zip" /Q
echo ---------------------------------------------------
echo Zipping up all files:
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "x64\Release\*.exe"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "x64\Release\*.dll"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "x64\Release\*.lib"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "*.dll"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "JMultiSliceLib.h"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "*.txt"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_BIN_Backup.zip" "*.pdf"
echo ---------------------------------------------------
echo Done.