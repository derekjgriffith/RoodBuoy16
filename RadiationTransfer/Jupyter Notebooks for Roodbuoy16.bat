REM This is a convenience batch file to run Jupyter notebooks in the root Anaconda environment.
REM Run this batch file (e.g. double-click from Windows file explorer) in order
REM to use the notebooks in the root Anaconda environment (not in a conda virtual environment) 
REM Deactivate any applicable virtual environment
REM deactivate is a batch script and must therefore be wrapped in a call command
call deactivate
jupyter kernelspec install-self --user
jupyter notebook



