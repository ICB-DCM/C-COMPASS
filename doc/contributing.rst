Contributing
============

We welcome contributions to C-COMPASS and encourage the community to participate in its development. Whether you are fixing bugs, adding new features, or improving documentation, your help is greatly appreciated.

To contribute, please follow these steps:

1. Fork the repository on GitHub and create a new branch for your work.
2. Make sure your changes adhere to the coding standards and include relevant tests, where applicable.
3. Submit a pull request with a clear description of the changes and the motivation behind them.
4. Ensure that your pull request is linked to any relevant issues or discussions.

Before starting major changes, it's a good idea to open an issue to discuss the proposed feature or bug fix. This helps avoid duplicate work and ensures your contributions are aligned with the project's goals. For additional guidance, please refer to our coding guidelines and the issue tracker on GitHub.

We appreciate your time and effort in making C-COMPASS even better!

Pre-commit Hooks
----------------

We use `pre-commit <https://github.com/pre-commit/pre-commit>`__ hooks to
ensure code quality and consistency. Pre-commit hooks automatically run checks
and formatting tools before each commit, helping to catch issues early.

To set up the pre-commit hooks in your local environment, follow these steps:

1. Install `pre-commit` if you haven't already:

   .. code-block:: sh

      pip install pre-commit
      # or install it together with the other development dependencies via
      pip install -e .[dev]

2. Navigate to the project directory and run:

   .. code-block:: sh

      pre-commit install

3. You're all set! The pre-commit hooks will now run automatically before each
   commit.

Building the Documentation
--------------------------

To build the documentation locally, you can use the following commands:

.. code-block:: sh

   # install tox if you haven't already
   pip install tox
   # build the documentation
   tox -e doc
   # this will generate the HTML documentation in doc/_build/html
   # you can open the documentation in your browser manually, or with
   python -c "import webbrowser; webbrowser.open('doc/_build/html/index.html')"

Release Process
---------------

To ensure a smooth release process, we follow these steps:

1. Update the release notes in ``doc/CHANGELOG.rst`` with the latest changes
   and version number.

2. Create a `new release <https://github.com/ICB-DCM/C-COMPASS/releases/new>`__
   via GitHub, following the Python versioning specifier scheme
   for the tag name (e.g., ``v1.0.0``). (For more information, see the
   `PyPA guide <https://packaging.python.org/en/latest/specifications/version-specifiers/#version-specifiers>`__.)

3. The release will be automatically built and published to PyPI and Zenodo
   by GitHub Actions. The documentation will also be updated on Read the Docs.
