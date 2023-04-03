# pybbat

## Initial Setup
This project is a template for a future app called `pybbat`.
Some modifications need to be made for it to do something useful!

-   **All sources & resources go into `pybbat/`**  
    Anything outside of this subdirectory does not get installed when the application is released.
-   **The `main` method in `pybbat/__main__.py` is the entrypoint**  
    You may modify this method to do whatever you'd like, however it must remain in `__main__.py`.
-   **Modify this readme to suit your application!**

## Applications as a Package
This application is a Python package, and needs to be treated like one during development. This means:
- Using Relative Imports &mdash; [Relative Import Info](https://realpython.com/absolute-vs-relative-python-imports/#relative-imports)
- Running the application as a module &mdash; `python -m pybbat`

None of the above points have any affect on the application after release; the application will be run as `pybbat` from the terminal.

## A Note on Virtual Environments
To ensure a reproducible build of this application, a virtual environment has been set up in the `venv/` directory. To set the packages in the environment, edit the appropriate `pybbat/requirements/*.in` setup file, then switch to that setup. For example, to set the production packages:
- Edit requirements/production.in
- Activate the environment: `source venv/bin/activate[.csh|.fish]`
- Compile the setup and synchronize your environment with it: `cadpip switch p`

## Creating a Release
1. Ensure all changes are commited using `git add` and `git commit`
2. Run `git release` from the terminal and follow prompts
3. Wait about a minute for your virtual environment to be assembled on the release server, and then run your application! (You can check the progress on the Gitlab continuous integration page of your project.)
