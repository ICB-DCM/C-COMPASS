Installation
============

See `https://github.com/ICB-DCM/C-COMPASS?tab=readme-ov-file#installation <https://github.com/ICB-DCM/C-COMPASS?tab=readme-ov-file#installation>`__ for the latest installation instructions.

Troubleshooting
---------------

- **SmartScreen Warning**: If Windows blocks the application via SmartScreen,
  it is because the software is unsigned. Please consult your IT department to
  bypass this restriction if necessary.

- **Long Path Issues on Windows**: If your system encounters long path errors,
  you can enable long paths in your registry:

  - Navigate to ``HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\FileSystem``.
  - Change the value of ``LongPathsEnabled`` from ``0`` to ``1``.
