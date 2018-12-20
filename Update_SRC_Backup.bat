@echo off
echo ---------------------------------------------------
echo Updating JMultiSliceLib C++ Source Backups.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
echo Removing previous zip files.
del "JMSLib_SRC_Backup.zip" /Q
echo ---------------------------------------------------
echo Zipping up all files:
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "src\*.h"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "src\*.cpp"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "src\cu\*.cuh"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "src\cu\*.cu"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "*.rtf"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "*.txt"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "*.manifest"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "*.bat"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "*.md"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "LICENSE"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_SRC_Backup.zip" "res\*.*"
echo ---------------------------------------------------
echo Done.