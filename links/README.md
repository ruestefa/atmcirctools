# Links directory

The `links/` directory contains soft-links to the source code of related projects, e.g., a library that is co-developed with the given application. This allows the library source files to be conveniently opened and edited in the same editor window as the application source files, while keeping the two repositories cleanly separated.

## Why not submodules?

The links mimic the "repo-within-a-repo" aspect of git submodules. However, git submodules have several drawbacks:

- Submodules are not compatible with conda builds. While it may be possible to specify a submodule as a pip dependency (provided it is even possible to load the submodule during a conda build), this breaks down as soon as the submodule has conda dependencies of its own. Those dependencies would need to be specified as dependencies of the main project, which erodes the separation between the two projects.

- Submodules could still be used during development, while specifying the submodule package as a regular dependency during installation. However, this requires the versions of the package in the submodule setup and in the dependencies to be kept in sync manually, which risks version incompatibilities.

- Submodules contain the whole project, not just the source code, which may include dozens of config files. These will all be visible when the submodule folder is opened in the file system tree of the editor, blowing it up in proportion. A minor annoyance, sure, but it negates many of the benefits of working on both projects in the same editor window.

## Why not separate repos?

That's actually exactly what's happening: The two repos are clones separately and bear no direct relation other than one being a dependency of the other.

To work on both source codes in the same browser window, the root directory of the editor's file system tree could just be set to the common root folder of woth repos. However, even if both are located right next to each other, this would still display both projects in full, including all top-level config files etc. (see last point of previous section).

If your editor allows multiple root nodes to be opened in the file system tree, this makes the links approach obsolete, of course, but not all editors support this.

## Caveats

The links-based approach is of course not without it's flaws:

- The relative symlinks assume a certain directory structure. If you organize your local repos differently, the links may be broken.

- Changes to a linked project have to be committed etc. separately. This is no different to git submodules, though, which also have to be entered in order to be accessed directly with git.
