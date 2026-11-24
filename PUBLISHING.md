# Source publication

The existing repository and remote are retained. `input/` experimental workbooks,
import caches and `output/` files have been removed from the current Git index,
without deleting the local files. `.gitignore` keeps them local. Small reference
LUTs and capacities are published explicitly under `input/reference/`.

Important: old commits still contain experimental files. An ordinary push of
that history distributes those files. No history was rewritten. For a release
that excludes them, use the prepared source-only ZIP, initialize a fresh
repository from that package, and review its contents before choosing the
publication destination. Replacing an existing remote/history is a separate
explicit action, not part of this preparation.

The source-only package contains neither `.git` nor raw workbooks/caches. It
includes the current source files, public references, tests and license notices.
The package is prepared from the current working tree, including existing user
changes, rather than an old committed revision.

Run `runModelOCPTests` in a clean MATLAB session before a release. Those tests do
not replace validation of full experimental fits. Code is CC BY-SA 4.0 with
inherited attribution, and the specified numerical references are CC BY 4.0.
Do not invent a release DOI or archived identifier; add actual metadata after
publication. No remote changes or pushes were performed by this preparation.
