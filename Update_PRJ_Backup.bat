@echo off
echo ---------------------------------------------------
echo Updating JMultiSliceLib C++ Project File Backups.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
echo Removing previous zip files.
del "JMSLib_PRJ_Backup.zip" /Q
echo ---------------------------------------------------
echo Zipping up all files:
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.sln"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.suo"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.vcxproj"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.vcxproj.*"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.txt"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.rtf"
"C:\Program Files\7-Zip\7z.exe" a "JMSLib_PRJ_Backup.zip" "*.pdf"
echo ---------------------------------------------------
echo Done.