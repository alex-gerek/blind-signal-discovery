<#
.SYNOPSIS
    Creates the "BSD Workstation" Windows Desktop shortcut pointing to the local
    launch_workstation.bat inside this repository, using the repository-stored icon.
#>
param(
    [string]$ShortcutName = "BSD Workstation"
)

$ErrorActionPreference = "Stop"

$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$LauncherPath = Join-Path $ProjectRoot "launch_workstation.bat"
$IconPath = Join-Path $ProjectRoot "app\assets\bsd_workstation.ico"

if (-not (Test-Path $LauncherPath)) {
    throw "Launcher not found: $LauncherPath"
}
if (-not (Test-Path $IconPath)) {
    throw "Icon not found: $IconPath"
}

$DesktopDir = [Environment]::GetFolderPath("Desktop")
if (-not (Test-Path $DesktopDir)) {
    throw "Windows Desktop directory not found: $DesktopDir"
}

$ShortcutPath = Join-Path $DesktopDir "$ShortcutName.lnk"
$WshShell = New-Object -ComObject WScript.Shell
$Shortcut = $WshShell.CreateShortcut($ShortcutPath)
$Shortcut.TargetPath = $LauncherPath
$Shortcut.WorkingDirectory = $ProjectRoot
$Shortcut.IconLocation = "$IconPath,0"
$Shortcut.Description = "Blind Signal Discovery (BSD) Local Measurement Workstation"
$Shortcut.WindowStyle = 1
$Shortcut.Save()

Write-Output "SHORTCUT_CREATED: $ShortcutPath"
Write-Output "TARGET: $($Shortcut.TargetPath)"
Write-Output "WORKDIR: $($Shortcut.WorkingDirectory)"
Write-Output "ICON: $($Shortcut.IconLocation)"
