# Launcher for agentcrabr Shiny app
# Run from PowerShell: .\launch.ps1

$AppDir = if ($PSScriptRoot) { $PSScriptRoot } else { (Get-Location).Path }
Set-Location $AppDir
Rscript "$AppDir/run.R"
