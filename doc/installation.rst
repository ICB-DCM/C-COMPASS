I. Installation
==============================

System Requirements
--------------------

- 64-bit Windows Operating System


Running the Software
--------------------

- Download the ZIP file from the repository or release section.
- Extract the ZIP file to any location on your machine.
- Navigate to the extracted folder.
- Double-click `C-CMPS.bat` to start the application.
- The software will initialize the portable Python environment and launch the GUI (this may take a few minutes).


Command-Line Usage (optional)
-----------------------------

You can also run the software via the command line:

> python CCMPS.py


Trouble-Shooting
----------------

- **SmartScreen Warning**: If Windows blocks the application via SmartScreen, it is because the software is unsigned. Please consult your IT department to bypass this restriction if necessary.
- **Long Path Issues on Windows**: If your system encounters long path errors, you can enable long paths in your registry:
   - Navigate to `HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\FileSystem`.
   - Set the value of **LongPathsEnabled** to `1`.
