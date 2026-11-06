# Launcher for agentcrabr Shiny app
# Double-click or run: powershell -File launch.ps1

$Rscript = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$AppDir = "C:\Users\Nathalie\projects\agentcrabr"
Set-Location $AppDir
& $Rscript -e "shiny::runApp('$AppDir/app.R', port=3840, launch.browser=TRUE)"
