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
   ```sh
   pip install pre-commit
   ```

2. Navigate to the project directory and run:
   ```sh
   pre-commit install
   ```

3. You're all set! The pre-commit hooks will now run automatically before each
   commit.
